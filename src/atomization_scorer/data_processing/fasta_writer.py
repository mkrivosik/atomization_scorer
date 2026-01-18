"""
fasta_writer.py

Provides functionality for writing FASTA sequence files.

Functions
---------
write_fasta : Writes sequences into a FASTA file.
"""

# ---------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------
from pathlib import Path
from typing import Dict

# ---------------------------------------------------------------------
# FASTA Writer
# ---------------------------------------------------------------------

def write_fasta(sequences: Dict[str, str], output_path: Path) -> Path:
    """
    Write sequences into a FASTA file.

    Parameters
    ----------
    sequences : dict
        Dictionary mapping FASTA headers to sequence strings.
    output_path : Path
        Path where the FASTA file should be written.

    Returns
    -------
    Path
        Path to the written FASTA file.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w") as file:
        for header, sequence in sequences.items():
            file.write(f">{header}\n{sequence}\n")

    return output_path
