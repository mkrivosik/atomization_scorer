"""
Tests for the compute_coverage_score() function.
"""

from pathlib import Path

import pytest

from atomization_scorer import compute_coverage_score


# --------------------------------------------------------------------------------------
# Test: returns a fraction between 0.0 and 1.0
# --------------------------------------------------------------------------------------
def test_compute_coverage_score_valid(mini_fasta: Path, mini_geese: Path):
    """compute_coverage_score should return a fraction between 0.0 and 1.0 for valid inputs."""
    score = compute_coverage_score(genomes_file=mini_fasta, atomization_file=mini_geese)
    assert 0.0 <= score <= 1.0


# --------------------------------------------------------------------------------------
# Test: returns 0.0 if genome has zero length
# --------------------------------------------------------------------------------------
def test_compute_coverage_score_zero_genome(mini_geese: Path, tmp_path: Path):
    """compute_coverage_score should return 0.0 if the genome has zero length."""
    empty_genomes = tmp_path / "empty.fasta"
    empty_genomes.write_text("")

    score = compute_coverage_score(genomes_file=empty_genomes, atomization_file=mini_geese)
    assert score == 0.0


# --------------------------------------------------------------------------------------
# Test: raises FileNotFoundError if genomes_file is missing
# --------------------------------------------------------------------------------------
def test_compute_coverage_missing_genomes_file(mini_geese: Path, tmp_path: Path):
    """compute_coverage_score should raise FileNotFoundError if genomes_file is missing."""
    missing_genomes = tmp_path / "missing.fasta"

    with pytest.raises(FileNotFoundError):
        compute_coverage_score(genomes_file=missing_genomes, atomization_file=mini_geese)


# --------------------------------------------------------------------------------------
# Test: raises FileNotFoundError if atomization_file is missing
# --------------------------------------------------------------------------------------
def test_compute_coverage_missing_atomization_file(mini_fasta: Path, tmp_path: Path):
    """compute_coverage_score should raise FileNotFoundError if atomization_file is missing."""
    missing_geese = tmp_path / "missing.geese"

    with pytest.raises(FileNotFoundError):
        compute_coverage_score(genomes_file=mini_fasta, atomization_file=missing_geese)


def test_compute_coverage_score_half_open_lengths(tmp_path: Path):
    """compute_coverage_score should use half-open interval lengths."""
    genomes_file = tmp_path / "genomes.fa"
    genomes_file.write_text(">sequence1\nAAAAAAAAAA\n>sequence2\nCCCCCCCCCC\n")

    atomization_file = tmp_path / "atoms.geese"
    atomization_file.write_text(
        "#name\tatom_nr\tclass\tstrand\tstart\tend\n"
        "sequence1\t1\t1\t+\t0\t4\n"
        "sequence2\t2\t2\t+\t5\t10\n"
    )

    score = compute_coverage_score(genomes_file=genomes_file, atomization_file=atomization_file)
    assert score == 9 / 20
