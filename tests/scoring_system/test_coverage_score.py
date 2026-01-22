"""
Tests for the compute_coverage_score() function.
"""

from pathlib import Path

import pytest

from atomization_scorer import compute_coverage_score


# ---------------------------------------------------------------------
# Test: returns a fraction between 0.0 and 1.0
# ---------------------------------------------------------------------
def test_compute_coverage_score_valid(mini_fasta: Path, mini_geese: Path):
    """compute_coverage_score should return a fraction between 0.0 and 1.0 for valid inputs."""
    score = compute_coverage_score(genomes_file=mini_fasta, atomization_file=mini_geese)
    assert 0.0 <= score <= 1.0


# ---------------------------------------------------------------------
# Test: returns 0.0 if genome has zero length
# ---------------------------------------------------------------------
def test_compute_coverage_score_zero_genome(mini_geese: Path, tmp_path: Path):
    """compute_coverage_score should return 0.0 if the genome has zero length."""
    empty_genomes = tmp_path / "empty.fasta"
    empty_genomes.write_text("")

    score = compute_coverage_score(genomes_file=empty_genomes, atomization_file=mini_geese)
    assert score == 0.0


# -----------------------------------------------------------------------------
# Test: raises FileNotFoundError if genomes_file is missing
# -----------------------------------------------------------------------------
def test_compute_coverage_missing_genomes_file(mini_geese: Path, tmp_path: Path):
    """compute_coverage_score should raise FileNotFoundError if genomes_file is missing."""
    missing_genomes = tmp_path / "missing.fasta"

    with pytest.raises(FileNotFoundError):
        compute_coverage_score(genomes_file=missing_genomes, atomization_file=mini_geese)


# -----------------------------------------------------------------------------
# Test: raises FileNotFoundError if atomization_file is missing
# -----------------------------------------------------------------------------
def test_compute_coverage_missing_atomization_file(mini_fasta: Path, tmp_path: Path):
    """compute_coverage_score should raise FileNotFoundError if atomization_file is missing."""
    missing_geese = tmp_path / "missing.geese"

    with pytest.raises(FileNotFoundError):
        compute_coverage_score(genomes_file=mini_fasta, atomization_file=missing_geese)
