"""
Tests for the coverage_score.py module.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
from pathlib import Path

import pandas as pd
import pytest
from Bio import SeqIO

from atomization_scorer import compute_coverage_score
from atomization_scorer.data_processing import read_geese


# --------------------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------------------
def _expected_coverage_score(fasta: Path, geese: Path) -> float:
    total_length = sum(len(record.seq) for record in SeqIO.parse(fasta, "fasta"))
    if total_length == 0:
        return 0.0
    atoms = read_geese(geese_file=geese).copy()
    atoms["start"] = pd.to_numeric(atoms["start"], errors="coerce")
    atoms["end"] = pd.to_numeric(atoms["end"], errors="coerce")
    atoms = atoms.dropna(subset=["start", "end"])
    return (atoms["end"] - atoms["start"]).sum() / total_length


def _expected_coverage_score_per_class(fasta: Path, geese: Path) -> list[dict[str, int | float]]:
    total_length = sum(len(record.seq) for record in SeqIO.parse(fasta, "fasta"))
    if total_length == 0:
        return []
    atoms = read_geese(geese_file=geese).copy()
    atoms["start"] = pd.to_numeric(atoms["start"], errors="coerce")
    atoms["end"] = pd.to_numeric(atoms["end"], errors="coerce")
    atoms = atoms.dropna(subset=["start", "end"])
    atoms["length"] = atoms["end"] - atoms["start"]
    results = []
    for atom_class, group in atoms.groupby("class"):
        results.append({"Class": int(str(atom_class)), "Coverage": float(group["length"].sum() / total_length)})
    return sorted(results, key=lambda entry: entry["Class"])


# --------------------------------------------------------------------------------------
# Test: computes the expected coverage score for valid inputs
# --------------------------------------------------------------------------------------
def test_compute_coverage_score(mini_fasta: Path, mini_geese: Path):
    """compute_coverage_score should return the expected coverage score for valid inputs."""
    score = compute_coverage_score(genomes_file=mini_fasta, atomization_file=mini_geese)
    assert score == _expected_coverage_score(mini_fasta, mini_geese)


# --------------------------------------------------------------------------------------
# Test: per-class coverage returns correct structure and values
# --------------------------------------------------------------------------------------
def test_compute_coverage_score_per_class(mini_fasta: Path, mini_geese: Path):
    """compute_coverage_score with per_class=True should return per-class coverage fractions."""
    result = compute_coverage_score(genomes_file=mini_fasta, atomization_file=mini_geese, per_class=True)
    expected = _expected_coverage_score_per_class(fasta=mini_fasta, geese=mini_geese)

    assert isinstance(result, list)
    assert all(isinstance(entry["Class"], int) for entry in result)
    assert all(0.0 <= entry["Coverage"] <= 1.0 for entry in result)
    assert result == sorted(result, key=lambda entry: entry["Class"])
    assert isinstance(expected, list)
    assert result == expected


# --------------------------------------------------------------------------------------
# Test: per-class coverages sum to overall coverage
# --------------------------------------------------------------------------------------
def test_compute_coverage_score_per_class_sums_to_overall(mini_fasta: Path, mini_geese: Path):
    """Per-class coverages should sum to the overall coverage score."""
    overall = compute_coverage_score(genomes_file=mini_fasta, atomization_file=mini_geese)
    per_class = compute_coverage_score(genomes_file=mini_fasta, atomization_file=mini_geese, per_class=True)

    assert isinstance(overall, float)
    assert isinstance(per_class, list)
    assert sum(entry["Coverage"] for entry in per_class) == pytest.approx(overall)


# --------------------------------------------------------------------------------------
# Test: returns 0.0 / [] if genome has zero length
# --------------------------------------------------------------------------------------
def test_compute_coverage_score_zero_genome(mini_geese: Path, tmp_path: Path):
    """compute_coverage_score should return 0.0 if the genome has zero length."""
    empty_genomes = tmp_path / "empty.fasta"
    empty_genomes.write_text("")

    score = compute_coverage_score(genomes_file=empty_genomes, atomization_file=mini_geese)
    assert score == 0.0


def test_compute_coverage_score_zero_genome_per_class(mini_geese: Path, tmp_path: Path):
    """compute_coverage_score with per_class=True should return [] if the genome has zero length."""
    empty_genomes = tmp_path / "empty.fasta"
    empty_genomes.write_text("")

    result = compute_coverage_score(genomes_file=empty_genomes, atomization_file=mini_geese, per_class=True)
    assert result == []


# --------------------------------------------------------------------------------------
# Test: raises FileNotFoundError if an input file is missing
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("missing_side", "expected_message"),
    [
        ("genomes", "Genomes file not found"),
        ("atomization", "Atomization file not found"),
    ],
)
def test_compute_coverage_score_missing_input_file(
    mini_fasta: Path,
    mini_geese: Path,
    tmp_path: Path,
    missing_side: str,
    expected_message: str,
):
    """compute_coverage_score should raise FileNotFoundError if an input file is missing."""
    genomes_file = mini_fasta
    atomization_file = mini_geese

    if missing_side == "genomes":
        genomes_file = tmp_path / "missing.fasta"
    else:
        atomization_file = tmp_path / "missing.geese"

    with pytest.raises(FileNotFoundError, match=expected_message):
        compute_coverage_score(genomes_file=genomes_file, atomization_file=atomization_file)
