"""
fasta_reader.py

Provides functionality for loading FASTA genomes sequences.

Functions
---------
read_fasta : Loads a FASTA file into a dictionary.
"""

# ---------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

# ---------------------------------------------------------------------
# FASTA Reader
# ---------------------------------------------------------------------

def read_fasta(fasta_file: Path) -> dict[str, Seq]:
    """
    Load a FASTA file and return its records as a dictionary.

    Parameters
    ----------
    fasta_file : Path
        Path to input FASTA genomes file.

    Raises
    ------
    FileNotFoundError
        Raised if the FASTA file does not exist.

    Returns
    -------
    dict
        Dictionary mapping sequence ID to sequence.
    """
    if not fasta_file.is_file():
        raise FileNotFoundError(f"FASTA file not found: {fasta_file}")

    return {records.id: records.seq for records in SeqIO.parse(fasta_file, "fasta")}
