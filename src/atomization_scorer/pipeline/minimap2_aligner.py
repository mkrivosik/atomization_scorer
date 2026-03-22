"""
minimap2_aligner.py

Utility function for running Minimap2 alignment.

Functions
---------
align_with_minimap2 : Align query sequences to target sequences using Minimap2.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
from __future__ import annotations
import logging
import subprocess
from pathlib import Path

# --------------------------------------------------------------------------------------
# Logging
# --------------------------------------------------------------------------------------
logger = logging.getLogger(__name__)

# --------------------------------------------------------------------------------------
# Minimap2 Alignment
# --------------------------------------------------------------------------------------

def align_with_minimap2(
    target: Path,
    query: Path,
    output_path: Path,
    preset: str = "asm20",
    emit_cigar: bool = True,
    secondary_ratio: float = 0.1,
) -> Path:
    """
    Align genomes sequences to atoms sequences using Minimap2 and generate PAF output.

    Parameters
    ----------
    target : Path
        Path to the target FASTA file.
    query : Path
         Path to the query FASTA file.
    output_path : Path
        Path where the resulting PAF alignment file should be saved.
    preset : str, optional, default="asm20"
        Minimap2 preset passed to ``-x``.
    emit_cigar : bool, optional, default=True
        Whether to include the ``-c`` option and emit CIGAR strings in PAF output.
    secondary_ratio : float, optional, default=0.1
        Minimap2 secondary-to-primary score ratio passed to ``-p``.

    Raises
    ------
    FileNotFoundError
        Raised if the target or query FASTA file does not exist, or if the
        minimap2 executable is not available on PATH.
    subprocess.CalledProcessError
        Raised if minimap2 fails during execution.
    ValueError
        Raised if minimap2 completes successfully but does not produce a
        non-empty PAF output file.

    Returns
    -------
    Path
        Generated PAF file path.
    """
    if not target.is_file():
        raise FileNotFoundError(f"Target FASTA file not found: {target}")
    if not query.is_file():
        raise FileNotFoundError(f"Query FASTA file not found: {query}")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Run minimap2
    command = [
        "minimap2",
        "-x", preset,
    ]

    if emit_cigar:
        command.append("-c")

    command.extend([
        "-p", str(secondary_ratio),
        str(target),
        str(query),
    ])

    logger.info("Running Minimap2 alignment:\n%s", " ".join(command))
    with output_path.open("w") as paf_file:
        try:
            subprocess.run(command, check=True, stdout=paf_file)
        except FileNotFoundError as error:
            raise FileNotFoundError("minimap2 executable not found on PATH.") from error

    if not output_path.is_file() or output_path.stat().st_size == 0:
        raise ValueError("Minimap2 did not produce a non-empty PAF output file.")

    logger.info("PAF alignment saved to %s", output_path)
    return output_path
