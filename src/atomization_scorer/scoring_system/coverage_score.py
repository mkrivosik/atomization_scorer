"""
coverage_score.py

Computes the coverage-based atomization score.

Modules
-------
compute_coverage_score  : Computes fraction of the genomes covered by atoms.
"""

# ---------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from atomization_scorer.data_processing import read_geese

# ---------------------------------------------------------------------
# Coverage Score Function
# ---------------------------------------------------------------------

def compute_coverage_score(
    genomes_file: Path,
    atomization_file: Path
) -> float:
    """
    Computes fraction of the genomes covered by atoms.

    Parameters
    ----------
    genomes_file : Path
        Input genomes FASTA file containing the genome sequences.
    atomization_file : Path
        Input GEESE file containing the predicted atomization.

    Raises
    ------
    FileNotFoundError
        Raised if genomes_file or atomization_file do not exist.

    Returns
    -------
    float
        Fraction of the genomes covered by atoms (0.0 to 1.0).
    """
    if not genomes_file.is_file():
        raise FileNotFoundError(f"Genomes file not found: {genomes_file}")
    if not atomization_file.is_file():
        raise FileNotFoundError(f"Atomization file not found: {atomization_file}")

    total_genomes_length = 0
    for record in SeqIO.parse(genomes_file, "fasta"):
        total_genomes_length += len(record.seq)

    if total_genomes_length == 0:
        return 0.0

    atoms_df = read_geese(geese_file=atomization_file)

    atoms_df["start"] = pd.to_numeric(atoms_df["start"], errors="coerce")
    atoms_df["end"] = pd.to_numeric(atoms_df["end"], errors="coerce")
    atoms_df = atoms_df.dropna(subset=["start", "end"])

    atoms_df["length"] = atoms_df["end"] - atoms_df["start"] + 1
    atoms_covered_length = atoms_df["length"].sum()

    return atoms_covered_length / total_genomes_length
