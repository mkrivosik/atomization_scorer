"""
Tests for the compute_overall_score() function.
"""

import pytest
from pathlib import Path
from atomization_scorer.scoring_system import compute_overall_score


# ---------------------------------------------------------------------
# Test: returns a float between 0.0 and 1.0
# ---------------------------------------------------------------------
def test_compute_overall_score_basic(mini_fasta: Path, mini_geese: Path, output_dir: Path, monkeypatch):
    """compute_overall_score should return a float between 0.0 and 1.0 for valid inputs."""
    monkeypatch.setattr(
        "atomization_scorer.scoring_system.overall_score.compute_alignment_score",
        lambda **kwargs: 0.8
    )
    monkeypatch.setattr(
        "atomization_scorer.scoring_system.overall_score.compute_coverage_score",
        lambda **kwargs: 0.9
    )

    overall_score = compute_overall_score(
        genomes_file=mini_fasta,
        atomization_file=mini_geese,
        output_directory=output_dir
    )

    assert overall_score == (0.8 ** 0.7) * (0.9 ** 0.3)


# ---------------------------------------------------------------------
# Test: returns 0.0 if scores are zero
# ---------------------------------------------------------------------
def test_compute_overall_score_zero(mini_fasta: Path, mini_geese: Path, output_dir: Path, monkeypatch):
    """compute_overall_score should return 0.0 if both alignment and coverage scores are zero."""
    monkeypatch.setattr(
        "atomization_scorer.scoring_system.overall_score.compute_alignment_score",
        lambda **kwargs: 0.0
    )
    monkeypatch.setattr(
        "atomization_scorer.scoring_system.overall_score.compute_coverage_score",
        lambda **kwargs: 0.0
    )

    overall = compute_overall_score(
        genomes_file=mini_fasta,
        atomization_file=mini_geese,
        output_directory=output_dir
    )
    assert overall == 0.0


# ---------------------------------------------------------------------
# Test: returns 1.0 if scores are one
# ---------------------------------------------------------------------
def test_compute_overall_score_one(mini_fasta: Path, mini_geese: Path, output_dir: Path, monkeypatch):
    """compute_overall_score should return 1.0 if both alignment and coverage scores are one."""
    monkeypatch.setattr(
        "atomization_scorer.scoring_system.overall_score.compute_alignment_score",
        lambda **kwargs: 1.0
    )
    monkeypatch.setattr(
        "atomization_scorer.scoring_system.overall_score.compute_coverage_score",
        lambda **kwargs: 1.0
    )

    overall = compute_overall_score(
        genomes_file=mini_fasta,
        atomization_file=mini_geese,
        output_directory=output_dir
    )
    assert overall == 1.0


# ---------------------------------------------------------------------
# Test: output_directory is created if it does not exist
# ---------------------------------------------------------------------
def test_compute_overall_score_creates_output_dir(mini_fasta: Path, mini_geese: Path, tmp_path: Path, monkeypatch):
    """compute_overall_score should create the output_directory if it does not exist."""
    non_existent_dir = tmp_path / "new_output_dir"

    monkeypatch.setattr(
        "atomization_scorer.scoring_system.overall_score.compute_alignment_score",
        lambda **kwargs: 0.5
    )
    monkeypatch.setattr(
        "atomization_scorer.scoring_system.overall_score.compute_coverage_score",
        lambda **kwargs: 0.5
    )

    assert not non_existent_dir.exists()

    compute_overall_score(
        genomes_file=mini_fasta,
        atomization_file=mini_geese,
        output_directory=non_existent_dir
    )

    assert non_existent_dir.exists()
    assert non_existent_dir.is_dir()


# ---------------------------------------------------------------------
# Test: raises FileNotFoundError if genomes_file does not exist
# ---------------------------------------------------------------------
def test_compute_overall_score_file_not_found(mini_geese: Path, output_dir: Path):
    """compute_overall_score should raise FileNotFoundError if genomes_file does not exist."""
    fake_fasta = output_dir / "not_exist.fa"

    with pytest.raises(FileNotFoundError):
        compute_overall_score(
            genomes_file=fake_fasta,
            atomization_file=mini_geese,
            output_directory=output_dir
        )


# ---------------------------------------------------------------------
# Test: raises FileNotFoundError if atomization_file does not exist
# ---------------------------------------------------------------------
def test_compute_overall_score_atomization_file_not_found(mini_fasta: Path, output_dir: Path):
    """compute_overall_score should raise FileNotFoundError if atomization_file does not exist."""
    fake_geese = output_dir / "not_exist.geese"

    with pytest.raises(FileNotFoundError):
        compute_overall_score(
            genomes_file=mini_fasta,
            atomization_file=fake_geese,
            output_directory=output_dir
        )
