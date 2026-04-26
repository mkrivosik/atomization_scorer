"""
Tests for the dotter_runner.py module.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
import subprocess
from pathlib import Path

import pytest

from atomization_scorer.diagnostics import dotter_runner as dr


# --------------------------------------------------------------------------------------
# Test: _get_display returns platform-appropriate DISPLAY value
# --------------------------------------------------------------------------------------
def test_get_display_on_macos_with_display_set(monkeypatch: pytest.MonkeyPatch):
    """_get_display should return host.docker.internal:<num> on macOS."""
    monkeypatch.setattr(dr.platform, "system", lambda: "Darwin")
    monkeypatch.setenv("DISPLAY", ":1")
    assert dr._get_display() == "host.docker.internal:1"


def test_get_display_on_macos_without_display_set(monkeypatch: pytest.MonkeyPatch):
    """_get_display should default to :0 and return host.docker.internal:0 on macOS."""
    monkeypatch.setattr(dr.platform, "system", lambda: "Darwin")
    monkeypatch.delenv("DISPLAY", raising=False)
    assert dr._get_display() == "host.docker.internal:0"


def test_get_display_on_linux_with_display_set(monkeypatch: pytest.MonkeyPatch):
    """_get_display should return the raw DISPLAY value on Linux."""
    monkeypatch.setattr(dr.platform, "system", lambda: "Linux")
    monkeypatch.setenv("DISPLAY", ":0")
    assert dr._get_display() == ":0"


def test_get_display_on_linux_without_display_set(monkeypatch: pytest.MonkeyPatch):
    """_get_display should return None when DISPLAY is unset on Linux."""
    monkeypatch.setattr(dr.platform, "system", lambda: "Linux")
    monkeypatch.delenv("DISPLAY", raising=False)
    assert dr._get_display() is None


# --------------------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------------------
def _create_anchor_directory(base_directory: Path, name: str = "rep_A") -> Path:
    """Create one anchor directory with minimal X/Y FASTA inputs."""
    anchor_directory = base_directory / name
    anchor_directory.mkdir(parents=True, exist_ok=True)
    (anchor_directory / "X.fasta").write_text(">anchor\nAAAA\n")
    (anchor_directory / "Y.fasta").write_text(">partner\nCCCC\n")
    return anchor_directory


# --------------------------------------------------------------------------------------
# Test: host paths are translated into container-visible paths
# --------------------------------------------------------------------------------------
def test_to_container_path_maps_files_under_mount_root(tmp_path: Path):
    """_to_container_path should map host files under the mount root into /work."""
    anchor_directory = _create_anchor_directory(tmp_path)

    assert dr._to_container_path(anchor_directory / "X.fasta", anchor_directory) == Path("/work/X.fasta")
    assert dr._to_container_path(anchor_directory / "dotter.pdf", anchor_directory) == Path("/work/dotter.pdf")


# --------------------------------------------------------------------------------------
# Test: nested host paths preserve their relative layout in the container
# --------------------------------------------------------------------------------------
def test_to_container_path_maps_nested_files_under_mount_root(tmp_path: Path):
    """_to_container_path should preserve nested relative paths under the mounted anchor directory."""
    anchor_directory = _create_anchor_directory(tmp_path / "anchors")
    nested_file = anchor_directory / "plots" / "dotter.pdf"
    nested_file.parent.mkdir(parents=True, exist_ok=True)
    nested_file.write_text("plot")

    assert dr._to_container_path(nested_file, anchor_directory) == Path("/work/plots/dotter.pdf")


# --------------------------------------------------------------------------------------
# Test: files outside the mount root are rejected
# --------------------------------------------------------------------------------------
def test_to_container_path_rejects_paths_outside_mount_root(tmp_path: Path):
    """_to_container_path should raise ValueError for paths outside the Docker mount root."""
    anchor_directory = _create_anchor_directory(tmp_path / "anchors")
    external_file = tmp_path / "external" / "X.fasta"
    external_file.parent.mkdir(parents=True, exist_ok=True)
    external_file.write_text(">x\nAAAA\n")

    with pytest.raises(ValueError, match="outside the Docker mount root"):
        dr._to_container_path(external_file, anchor_directory)


# --------------------------------------------------------------------------------------
# Test: unsupported output format is rejected
# --------------------------------------------------------------------------------------
def test_run_dotter_for_anchor_unsupported_output_format(tmp_path: Path):
    """run_dotter_for_anchor should raise ValueError for unsupported output formats."""
    anchor_directory = _create_anchor_directory(tmp_path)

    with pytest.raises(ValueError, match="Unsupported Dotter output format"):
        dr.run_dotter_for_anchor(
            anchor_directory=anchor_directory,
            output_format="jpg",
        )


# --------------------------------------------------------------------------------------
# Test: missing anchor inputs are rejected
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("missing_side", "expected_message"),
    [
        ("anchor_directory", "Anchor directory not found"),
        ("x_fasta", "Anchor X FASTA file not found"),
        ("y_fasta", "Anchor Y FASTA file not found"),
    ],
)
def test_run_dotter_for_anchor_missing_input(
    tmp_path: Path,
    missing_side: str,
    expected_message: str,
):
    """run_dotter_for_anchor should raise FileNotFoundError when required inputs are missing."""
    anchor_directory = tmp_path / "rep_A"
    if missing_side != "anchor_directory":
        anchor_directory.mkdir()
    if missing_side not in {"anchor_directory", "x_fasta"}:
        (anchor_directory / "X.fasta").write_text(">anchor\nAAAA\n")
    if missing_side not in {"anchor_directory", "y_fasta"}:
        (anchor_directory / "Y.fasta").write_text(">partner\nCCCC\n")

    with pytest.raises(FileNotFoundError, match=expected_message):
        dr.run_dotter_for_anchor(anchor_directory=anchor_directory)


# --------------------------------------------------------------------------------------
# Test: runs Dotter for one anchor with supported output formats on Linux
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize("output_format", ["png", "svg", "pdf"])
def test_run_dotter_for_anchor_on_linux(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    output_format: str,
):
    """On Linux without DISPLAY or X11 socket, run_dotter_for_anchor should build the minimal Docker command."""
    anchor_directory = _create_anchor_directory(tmp_path)
    calls = []

    def fake_run(*args, **kwargs):
        calls.append((args[0], kwargs))
        return subprocess.CompletedProcess(args[0], 0, stdout="dotter-output\n", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)
    monkeypatch.setattr(dr.platform, "system", lambda: "Linux")
    monkeypatch.delenv("DISPLAY", raising=False)
    monkeypatch.setattr(dr, "X11_SOCKET_DIRECTORY", tmp_path / "missing-x11")

    dr.run_dotter_for_anchor(
        anchor_directory=anchor_directory,
        extra_args=["-v"],
        output_format=output_format,
    )

    assert calls == [
        (
            [
                "docker",
                "run",
                "--rm",
                "-v",
                f"{anchor_directory.resolve()}:/work",
                dr.DOTTER_IMAGE,
                "dotter",
                "-e",
                f"/work/dotter.{output_format}",
                "-v",
                "/work/Y.fasta",
                "/work/X.fasta",
            ],
            {"check": True, "capture_output": True, "text": True},
        )
    ]


# --------------------------------------------------------------------------------------
# Test: on macOS DISPLAY is always injected even when unset, socket is never mounted
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize("output_format", ["png", "svg", "pdf"])
def test_run_dotter_for_anchor_on_macos(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    output_format: str,
):
    """On macOS without DISPLAY or X11 socket, run_dotter_for_anchor should always inject
    DISPLAY=host.docker.internal:0."""
    anchor_directory = _create_anchor_directory(tmp_path)
    calls = []

    def fake_run(*args, **kwargs):
        calls.append((args[0], kwargs))
        return subprocess.CompletedProcess(args[0], 0, stdout="dotter-output\n", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)
    monkeypatch.setattr(dr.platform, "system", lambda: "Darwin")
    monkeypatch.delenv("DISPLAY", raising=False)
    monkeypatch.setattr(dr, "X11_SOCKET_DIRECTORY", tmp_path / "missing-x11")

    dr.run_dotter_for_anchor(
        anchor_directory=anchor_directory,
        extra_args=["-v"],
        output_format=output_format,
    )

    assert calls == [
        (
            [
                "docker",
                "run",
                "--rm",
                "--env",
                "DISPLAY=host.docker.internal:0",
                "-v",
                f"{anchor_directory.resolve()}:/work",
                dr.DOTTER_IMAGE,
                "dotter",
                "-e",
                f"/work/dotter.{output_format}",
                "-v",
                "/work/Y.fasta",
                "/work/X.fasta",
            ],
            {"check": True, "capture_output": True, "text": True},
        )
    ]


# --------------------------------------------------------------------------------------
# Test: X11 settings are forwarded on Linux (socket mounted, DISPLAY passed as-is)
# --------------------------------------------------------------------------------------
def test_run_dotter_for_anchor_forwards_x11_runtime_on_linux(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    """On Linux, run_dotter_for_anchor should forward DISPLAY as-is and mount the X11 socket."""
    anchor_directory = _create_anchor_directory(tmp_path)
    x11_socket_directory = tmp_path / ".X11-unix"
    x11_socket_directory.mkdir()
    calls = []

    def fake_run(*args, **kwargs):
        calls.append((args[0], kwargs))
        return subprocess.CompletedProcess(args[0], 0, stdout="dotter-output\n", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)
    monkeypatch.setattr(dr.platform, "system", lambda: "Linux")
    monkeypatch.setenv("DISPLAY", ":0")
    monkeypatch.setattr(dr, "X11_SOCKET_DIRECTORY", x11_socket_directory)

    dr.run_dotter_for_anchor(anchor_directory=anchor_directory)

    assert calls == [
        (
            [
                "docker",
                "run",
                "--rm",
                "--env",
                "DISPLAY=:0",
                "-v",
                f"{x11_socket_directory}:{x11_socket_directory}",
                "-v",
                f"{anchor_directory.resolve()}:/work",
                dr.DOTTER_IMAGE,
                "dotter",
                "-e",
                "/work/dotter.pdf",
                "/work/Y.fasta",
                "/work/X.fasta",
            ],
            {"check": True, "capture_output": True, "text": True},
        )
    ]


# --------------------------------------------------------------------------------------
# Test: X11 settings on macOS use host.docker.internal and skip socket mount
# --------------------------------------------------------------------------------------
def test_run_dotter_for_anchor_forwards_x11_runtime_on_macos(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    """On macOS, run_dotter_for_anchor should route DISPLAY via host.docker.internal and not mount the socket."""
    anchor_directory = _create_anchor_directory(tmp_path)
    x11_socket_directory = tmp_path / ".X11-unix"
    x11_socket_directory.mkdir()
    calls = []

    def fake_run(*args, **kwargs):
        calls.append((args[0], kwargs))
        return subprocess.CompletedProcess(args[0], 0, stdout="dotter-output\n", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)
    monkeypatch.setattr(dr.platform, "system", lambda: "Darwin")
    monkeypatch.setenv("DISPLAY", ":0")
    monkeypatch.setattr(dr, "X11_SOCKET_DIRECTORY", x11_socket_directory)

    dr.run_dotter_for_anchor(anchor_directory=anchor_directory)

    assert calls == [
        (
            [
                "docker",
                "run",
                "--rm",
                "--env",
                "DISPLAY=host.docker.internal:0",
                "-v",
                f"{anchor_directory.resolve()}:/work",
                dr.DOTTER_IMAGE,
                "dotter",
                "-e",
                "/work/dotter.pdf",
                "/work/Y.fasta",
                "/work/X.fasta",
            ],
            {"check": True, "capture_output": True, "text": True},
        )
    ]


# --------------------------------------------------------------------------------------
# Test: Dotter runtime failures propagate to the caller
# --------------------------------------------------------------------------------------
def test_run_dotter_for_anchor_propagates_called_process_error(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    """run_dotter_for_anchor should propagate a subprocess.CalledProcessError on Docker failure."""
    anchor_directory = _create_anchor_directory(tmp_path)

    def fake_run(*_args, **_kwargs):
        raise subprocess.CalledProcessError(1, ["docker", "run"], stderr="dotter failed")

    monkeypatch.setattr(subprocess, "run", fake_run)

    with pytest.raises(subprocess.CalledProcessError):
        dr.run_dotter_for_anchor(anchor_directory=anchor_directory)


# --------------------------------------------------------------------------------------
# Test: missing Dotter executable is rejected
# --------------------------------------------------------------------------------------
def test_run_dotter_for_anchor_missing_dotter_executable(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    """run_dotter_for_anchor should raise FileNotFoundError if Docker is not on PATH."""
    anchor_directory = _create_anchor_directory(tmp_path)

    def fake_run(*_args, **_kwargs):
        raise FileNotFoundError("docker")

    monkeypatch.setattr(subprocess, "run", fake_run)

    with pytest.raises(FileNotFoundError, match="docker executable not found on PATH"):
        dr.run_dotter_for_anchor(anchor_directory=anchor_directory)


# --------------------------------------------------------------------------------------
# Test: missing anchors directory is rejected
# --------------------------------------------------------------------------------------
def test_run_dotter_for_anchors_missing_directory(tmp_path: Path):
    """run_dotter_for_anchors should raise FileNotFoundError if the anchors directory does not exist."""
    with pytest.raises(FileNotFoundError, match="Anchors directory not found"):
        dr.run_dotter_for_anchors(anchors_directory=tmp_path / "missing")


# --------------------------------------------------------------------------------------
# Test: runs Dotter across all valid anchor directories
# --------------------------------------------------------------------------------------
def test_run_dotter_for_anchors(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    """run_dotter_for_anchors should process all anchor directories containing X/Y FASTA files."""
    anchors_directory = tmp_path / "anchors"
    _create_anchor_directory(anchors_directory, "rep_A")
    _create_anchor_directory(anchors_directory, "rep_B")
    skipped_directory = anchors_directory / "rep_C"
    skipped_directory.mkdir(parents=True, exist_ok=True)
    (skipped_directory / "X.fasta").write_text(">anchor\nAAAA\n")

    calls = []

    def fake_run(*args, **kwargs):
        calls.append((args[0], kwargs))
        return subprocess.CompletedProcess(args[0], 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)
    monkeypatch.setattr(dr.platform, "system", lambda: "Linux")
    monkeypatch.delenv("DISPLAY", raising=False)
    monkeypatch.setattr(dr, "X11_SOCKET_DIRECTORY", tmp_path / "missing-x11")

    dr.run_dotter_for_anchors(
        anchors_directory=anchors_directory,
    )

    assert len(calls) == 2
    assert calls[0][0][6:9] == ["dotter", "-e", "/work/dotter001.pdf"]
    assert calls[1][0][6:9] == ["dotter", "-e", "/work/dotter002.pdf"]


# --------------------------------------------------------------------------------------
# Test: multi-anchor runner forwards extra arguments and output format
# --------------------------------------------------------------------------------------
def test_run_dotter_for_anchors_forwards_extra_args_and_output_format(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    """run_dotter_for_anchors should pass extra_args and output_format to each anchor run."""
    anchors_directory = tmp_path / "anchors"
    anchor_a = _create_anchor_directory(anchors_directory, "rep_A")
    anchor_b = _create_anchor_directory(anchors_directory, "rep_B")
    calls = []

    def fake_run_dotter_for_anchor(
        anchor_directory: Path,
        extra_args=None,
        output_format: str = "pdf",
        output_stem: str = "dotter",
    ):
        calls.append((anchor_directory, extra_args, output_format, output_stem))

    monkeypatch.setattr(dr, "run_dotter_for_anchor", fake_run_dotter_for_anchor)

    dr.run_dotter_for_anchors(
        anchors_directory=anchors_directory,
        extra_args=["-v", "-q"],
        output_format="svg",
    )

    assert calls == [
        (anchor_a, ["-v", "-q"], "svg", "dotter001"),
        (anchor_b, ["-v", "-q"], "svg", "dotter002"),
    ]
