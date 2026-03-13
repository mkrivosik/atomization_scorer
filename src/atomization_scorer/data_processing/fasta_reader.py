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
from Bio.SeqIO.FastaIO import SimpleFastaParser

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
    ValueError
        Raised if the FASTA file contains duplicate IDs or malformed FASTA content.

    Returns
    -------
    dict
        Dictionary mapping sequence ID to sequence.
    """
    if not fasta_file.is_file():
        raise FileNotFoundError(f"FASTA file not found: {fasta_file}")

    try:
        with fasta_file.open() as file:
            parsed_headers = [title for title, _ in SimpleFastaParser(file)]
    except ValueError as error:
        raise ValueError(f"Malformed FASTA file: {fasta_file}") from error

    if any(not header for header in parsed_headers):
        raise ValueError(f"Malformed FASTA file: {fasta_file}")

    records = {}
    seen_ids = set()
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id in seen_ids:
            raise ValueError(f"Duplicate FASTA ID found: {record.id}")
        seen_ids.add(record.id)
        records[record.id] = record.seq

    return records
