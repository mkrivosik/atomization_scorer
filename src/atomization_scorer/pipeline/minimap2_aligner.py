"""
minimap2_aligner.py

Utility function for running Minimap2 alignment.

Functions
---------
align_with_minimap2 : Align query sequences to target sequences using Minimap2.
"""

# ---------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------
import subprocess
from pathlib import Path

# ---------------------------------------------------------------------
# Minimap2 Alignment
# ---------------------------------------------------------------------

def align_with_minimap2(
    target: Path,
    query: Path,
    output_path: Path,
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
        Path where resulting PAF alignment file should be saved.

    Raises
    ------
    FileNotFoundError
        Raised if target or query file do not exist.
    subprocess.CalledProcessError
        Raised if minimap2 fails during execution.

    Returns
    -------
    Path
        Path to generated PAF file.
    """
    if not target.is_file():
        raise FileNotFoundError(f"Target FASTA file not found: {target}")
    if not query.is_file():
        raise FileNotFoundError(f"Query FASTA file not found: {query}")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Run minimap2
    command = [
        "minimap2",
        "-x", "asm20",
        "-c",
        "-p", "0.1",
        str(target),
        str(query),
    ]

    print(f"Running Minimap2 alignment:\n{' '.join(command)}")
    with output_path.open("w") as paf_file:
        subprocess.run(command, check=True, stdout=paf_file)

    print(f"PAF alignment saved to {output_path}")
    return output_path
