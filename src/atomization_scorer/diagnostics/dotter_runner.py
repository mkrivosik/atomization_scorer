"""
dotter_runner.py

Provides utilities for running Dotter in Docker on generated
anchor-vs.-partners FASTA inputs.

Functions
---------
_to_container_path      : Map one host path into the container work directory.
_build_dotter_command   : Build the Docker command used to run Dotter for one anchor.
run_dotter_for_anchor   : Run Dotter for one anchor directory and export a plot.
run_dotter_for_anchors  : Run Dotter for all anchor directories and export plots.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
from __future__ import annotations
import logging
import os
import platform
import subprocess
from pathlib import Path

# --------------------------------------------------------------------------------------
# Logging
# --------------------------------------------------------------------------------------
log = logging.getLogger(__name__)

# --------------------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------------------
DOTTER_IMAGE                = "nathanhaigh/seqtools:4.44.1"
CONTAINER_WORKDIR           = Path("/work")
SUPPORTED_OUTPUT_FORMATS    = {"png", "svg", "pdf"}
X11_SOCKET_DIRECTORY        = Path("/tmp/.X11-unix")


# --------------------------------------------------------------------------------------
# Display Helpers
# --------------------------------------------------------------------------------------
def _get_display() -> str | None:
    """
    Return the DISPLAY value to pass into the Docker container.

    On macOS Docker runs inside a Linux VM and cannot reach the host Unix socket
    at /tmp/.X11-unix; the display is routed via host.docker.internal instead.
    On Linux the host DISPLAY value is returned as-is.

    Returns
    -------
    str or None
        On macOS, host.docker.internal:<display_num> derived from DISPLAY (defaults
        to :0 if unset). On Linux, the raw DISPLAY environment variable, or None
        if unset.
    """
    if platform.system() == "Darwin":
        display = os.environ.get("DISPLAY", ":0")
        display_num = display.lstrip(":")
        return f"host.docker.internal:{display_num}"
    return os.environ.get("DISPLAY")


# --------------------------------------------------------------------------------------
# Path Mapping Helpers
# --------------------------------------------------------------------------------------
def _to_container_path(host_path: Path, mount_root: Path) -> Path:
    """
    Map one host path under the Docker mount root into the container work directory.

    Parameters
    ----------
    host_path : Path
        Host path to map into the container.
    mount_root : Path
        Host directory mounted into the container at /work.

    Raises
    ------
    ValueError
        Raised if the host path is not located under the mount root.

    Returns
    -------
    Path
        Container-visible path under /work.
    """
    resolved_host_path = host_path.resolve()
    resolved_mount_root = mount_root.resolve()

    try:
        relative_path = resolved_host_path.relative_to(resolved_mount_root)
    except ValueError as error:
        raise ValueError(
            f"Path {resolved_host_path} is outside the Docker mount root {resolved_mount_root}."
        ) from error

    return CONTAINER_WORKDIR / relative_path


def _build_dotter_command(
    anchor_directory: Path,
    x_fasta: Path,
    y_fasta: Path,
    output_file: Path,
    extra_args: list[str] | None = None,
    image: str = DOTTER_IMAGE,
) -> list[str]:
    """
    Build the Docker command used to run Dotter for one anchor directory.

    Parameters
    ----------
    anchor_directory : Path
        Host directory mounted into the container.
    x_fasta : Path
        Host path to the X FASTA file under anchor_directory.
    y_fasta : Path
        Host path to the Y FASTA file under anchor_directory.
    output_file : Path
        Host path to the exported Dotter output under anchor_directory.
    extra_args : list[str] or None, optional, default=None
        Additional command-line arguments passed to Dotter.
    image : str, optional, default=DOTTER_IMAGE
        Docker image containing Dotter.

    Returns
    -------
    list[str]
        Docker command for batch Dotter execution.
    """
    resolved_anchor_directory = anchor_directory.resolve()
    container_x = _to_container_path(x_fasta, resolved_anchor_directory)
    container_y = _to_container_path(y_fasta, resolved_anchor_directory)
    container_output = _to_container_path(output_file, resolved_anchor_directory)

    command = [
        "docker",
        "run",
        "--rm",
    ]
    display = _get_display()
    if display:
        command.extend(["--env", f"DISPLAY={display}"])
    if platform.system() != "Darwin" and X11_SOCKET_DIRECTORY.is_dir():
        command.extend(["-v", f"{X11_SOCKET_DIRECTORY}:{X11_SOCKET_DIRECTORY}"])
    command.extend([
        "-v",
        f"{resolved_anchor_directory}:{CONTAINER_WORKDIR}",
        image,
        "dotter",
        "-e",
        str(container_output),
    ])
    if extra_args:
        command.extend(extra_args)
    command.extend([str(container_y), str(container_x)])
    return command


# --------------------------------------------------------------------------------------
# Single-Anchor Dotter Runner
# --------------------------------------------------------------------------------------
def run_dotter_for_anchor(
    anchor_directory: Path,
    extra_args: list[str] | None = None,
    output_format: str = "pdf",
    output_stem: str = "dotter",
) -> None:
    """
    Run Dotter for one anchor directory.

    Parameters
    ----------
    anchor_directory : Path
        Anchor directory containing X.fasta and Y.fasta.
    extra_args : list[str] or None, optional, default=None
        Additional command-line arguments passed to Dotter.
    output_format : str, optional, default="pdf"
        Output format written into the anchor directory.
    output_stem : str, optional, default="dotter"
        Basename of the exported Dotter output file without extension.

    Raises
    ------
    FileNotFoundError
        Raised if the anchor directory, X.fasta, or Y.fasta does not exist,
        or if the Docker executable is not available on PATH.
    ValueError
        Raised if the requested output format is not supported.
    subprocess.CalledProcessError
        Raised if Dotter fails during execution.

    Returns
    -------
    None
    """
    if not anchor_directory.is_dir():
        raise FileNotFoundError(f"Anchor directory not found: {anchor_directory}")

    x_fasta = anchor_directory / "X.fasta"
    y_fasta = anchor_directory / "Y.fasta"
    if not x_fasta.is_file():
        raise FileNotFoundError(f"Anchor X FASTA file not found: {x_fasta}")
    if not y_fasta.is_file():
        raise FileNotFoundError(f"Anchor Y FASTA file not found: {y_fasta}")
    if output_format not in SUPPORTED_OUTPUT_FORMATS:
        raise ValueError(
            f"Unsupported Dotter output format: {output_format}. "
            f"Supported formats: {sorted(SUPPORTED_OUTPUT_FORMATS)}"
        )

    output_file = anchor_directory / f"{output_stem}.{output_format}"
    command = _build_dotter_command(
        anchor_directory=anchor_directory,
        x_fasta=x_fasta,
        y_fasta=y_fasta,
        output_file=output_file,
        extra_args=extra_args,
    )
    log.info("Running Dotter for anchor directory %s:\n%s", anchor_directory, " ".join(command))

    try:
        subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError as error:
        raise FileNotFoundError("docker executable not found on PATH.") from error

    return None


# --------------------------------------------------------------------------------------
# Multi-Anchor Dotter Runner
# --------------------------------------------------------------------------------------
def run_dotter_for_anchors(
    anchors_directory: Path,
    extra_args: list[str] | None = None,
    output_format: str = "pdf",
) -> None:
    """
    Run Dotter for all anchor directories.

    Parameters
    ----------
    anchors_directory : Path
        Directory containing per-anchor subdirectories with X.fasta and Y.fasta.
    extra_args : list[str] or None, optional, default=None
        Additional command-line arguments passed to Dotter.
    output_format : str, optional, default="pdf"
        Output format written into each anchor directory.

    Raises
    ------
    FileNotFoundError
        Raised if the anchors directory does not exist.
    ValueError
        Raised if the requested output format is not supported.

    Returns
    -------
    None
    """
    if not anchors_directory.is_dir():
        raise FileNotFoundError(f"Anchors directory not found: {anchors_directory}")

    dotter_index = 1
    for anchor_directory in sorted(path for path in anchors_directory.iterdir() if path.is_dir()):
        if not (anchor_directory / "X.fasta").is_file():
            continue
        if not (anchor_directory / "Y.fasta").is_file():
            continue

        run_dotter_for_anchor(
            anchor_directory=anchor_directory,
            extra_args=extra_args,
            output_format=output_format,
            output_stem=f"dotter{dotter_index:03d}",
        )
        dotter_index += 1

    return None
