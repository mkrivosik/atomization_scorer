"""
coverage_score.py

Computes the coverage-based atomization score.

Functions
---------
compute_coverage_score  : Computes a fraction of the genomes covered by atoms.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
from __future__ import annotations
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from atomization_scorer.data_processing import read_geese


# --------------------------------------------------------------------------------------
# Coverage Score Function
# --------------------------------------------------------------------------------------
def compute_coverage_score(
    genomes_file: Path,
    atomization_file: Path,
    per_class: bool = False,
) -> float | list[dict[str, int | float]]:
    """
    Computes a fraction of the genomes covered by atoms.

    Parameters
    ----------
    genomes_file : Path
        Input genomes FASTA file containing the genome sequences.
    atomization_file : Path
        Input GEESE file containing the predicted atomization.
    per_class : bool, optional, default=False
        If True, returns per-class coverage fractions instead of the overall fraction.

    Raises
    ------
    FileNotFoundError
        Raised if genomes_file or atomization_file do not exist.

    Returns
    -------
    float or list[dict[str, int | float]]
        If per_class is False, returns the overall fraction of genomes covered (0.0 to 1.0).
        If per_class is True, returns a list of dictionaries sorted by class, each containing:
            "Class": int -> atomization class,
            "Coverage": float -> fraction of total genome length covered by that class (0.0 to 1.0).
    """
    if not genomes_file.is_file():
        raise FileNotFoundError(f"Genomes file not found: {genomes_file}")
    if not atomization_file.is_file():
        raise FileNotFoundError(f"Atomization file not found: {atomization_file}")

    total_genomes_length = 0
    for record in SeqIO.parse(genomes_file, "fasta"):
        total_genomes_length += len(record.seq)

    if total_genomes_length == 0:
        return [] if per_class else 0.0

    atoms_df = read_geese(geese_file=atomization_file).copy()

    atoms_df["start"] = pd.to_numeric(atoms_df["start"], errors="coerce")
    atoms_df["end"] = pd.to_numeric(atoms_df["end"], errors="coerce")
    atoms_df = atoms_df.dropna(subset=["start", "end"]).copy()

    # Based on the IMP atom definition, atoms are assumed not to overlap.
    atoms_df["length"] = atoms_df["end"] - atoms_df["start"]

    if per_class:
        results = []
        for atom_class, group in atoms_df.groupby("class"):
            coverage = float(group["length"].sum() / total_genomes_length)
            results.append({"Class": int(str(atom_class)), "Coverage": coverage})
        return sorted(results, key=lambda entry: entry["Class"])

    return float(atoms_df["length"].sum() / total_genomes_length)
