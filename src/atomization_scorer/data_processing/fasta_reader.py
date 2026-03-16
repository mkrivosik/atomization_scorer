"""
fasta_reader.py

Provides functionality for loading FASTA genomes sequences.

Functions
---------
read_fasta : Loads a FASTA file into a dictionary.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser

# --------------------------------------------------------------------------------------
# FASTA Reader
# --------------------------------------------------------------------------------------

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

    records = {}
    try:
        with fasta_file.open() as file:
            for title, sequence in SimpleFastaParser(file):
                if not title:
                    raise ValueError(f"Malformed FASTA file: {fasta_file}")

                record_id = title.split(None, 1)[0]
                if record_id in records:
                    raise ValueError(f"Duplicate FASTA ID found: {record_id}")

                records[record_id] = Seq(sequence)
    except ValueError as error:
        if str(error).startswith("Duplicate FASTA ID found:"):
            raise
        if str(error).startswith("Malformed FASTA file:"):
            raise
        raise ValueError(f"Malformed FASTA file: {fasta_file}") from error

    return records
