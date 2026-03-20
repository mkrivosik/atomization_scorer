"""
Tests for the compute_coverage_score() function.
"""

from pathlib import Path

import pytest

from atomization_scorer import compute_coverage_score


# --------------------------------------------------------------------------------------
# Test: computes the expected coverage score for valid inputs
# --------------------------------------------------------------------------------------
def test_compute_coverage_score(mini_fasta: Path, mini_geese: Path):
    """compute_coverage_score should return the expected coverage score for valid inputs."""
    score = compute_coverage_score(genomes_file=mini_fasta, atomization_file=mini_geese)
    assert score == 4981779 / 5004626


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
