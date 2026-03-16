"""
fasta_writer.py

Provides functionality for writing FASTA sequence files.

Functions
---------
write_fasta : Writes sequences into a FASTA file.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
from pathlib import Path

from Bio.Seq import Seq

# --------------------------------------------------------------------------------------
# FASTA Writer
# --------------------------------------------------------------------------------------

def write_fasta(sequences: dict[str, Seq], output_path: Path) -> Path:
    """
    Write sequences into a FASTA file.

    Parameters
    ----------
    sequences : dict
        Dictionary mapping FASTA headers to sequence values.
    output_path : Path
        Path where the FASTA file should be written.

    Raises
    ------
    ValueError
        Raised if a FASTA header contains unsafe characters or a sequence is empty.

    Returns
    -------
    Path
        Written FASTA file path.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w") as file:
        for header, sequence in sequences.items():
            if not header:
                raise ValueError("FASTA header must not be empty.")
            if ">" in header:
                raise ValueError(f"FASTA header contains '>': {header}")
            if "\n" in header or "\r" in header:
                raise ValueError(f"FASTA header contains newline characters: {header}")

            sequence_text = str(sequence)
            if not sequence_text:
                raise ValueError(f"FASTA sequence must not be empty for header: {header}")

            file.write(f">{header}\n")
            for start in range(0, len(sequence_text), 80):
                file.write(f"{sequence_text[start:start + 80]}\n")

    return output_path
