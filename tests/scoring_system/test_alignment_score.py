"""
Tests for compute_alignment_score() function.
"""

from pathlib import Path

import pytest

from atomization_scorer.scoring_system import compute_alignment_score


# -----------------------------------------------------------------------------
# Test: basic base-level score computation
# -----------------------------------------------------------------------------
def test_compute_alignment_score_base(mini_fasta: Path, mini_geese: Path, output_dir: Path, monkeypatch):
    """compute_alignment_score should return correct base-level score."""

    monkeypatch.setattr(
        "atomization_scorer.scoring_system.alignment_score.compute_true_alignment",
        lambda *args, **kwargs: mini_geese
    )

    monkeypatch.setattr(
        "atomization_scorer.scoring_system.alignment_score.compute_base_level_metrics",
        lambda *args, **kwargs: 0.85
    )

    score = compute_alignment_score(
        genomes_file=mini_fasta,
        atomization_file=mini_geese,
        output_directory=output_dir,
        level="base",
        per_class=False
    )

    assert isinstance(score, float)
    assert 0.0 <= score <= 1.0
    assert score == 0.85


# -----------------------------------------------------------------------------
# Test: interval-level score computation
# -----------------------------------------------------------------------------
def test_compute_alignment_score_interval(mini_fasta: Path, mini_geese: Path, output_dir: Path, monkeypatch):
    """compute_alignment_score should return correct interval-level score."""

    monkeypatch.setattr(
        "atomization_scorer.scoring_system.alignment_score.compute_true_alignment",
        lambda *args, **kwargs: mini_geese
    )

    monkeypatch.setattr(
        "atomization_scorer.scoring_system.alignment_score.compute_interval_level_metrics",
        lambda *args, **kwargs: 0.75
    )

    score = compute_alignment_score(
        genomes_file=mini_fasta,
        atomization_file=mini_geese,
        output_directory=output_dir,
        level="interval",
        per_class=True,
        min_overlap_ratio=0.9
    )

    assert isinstance(score, float)
    assert 0.0 <= score <= 1.0
    assert score == 0.75


# -----------------------------------------------------------------------------
# Test: raises FileNotFoundError if genomes_file is missing
# -----------------------------------------------------------------------------
def test_compute_alignment_score_missing_genomes(mini_geese: Path, tmp_path: Path, output_dir: Path):
    """compute_alignment_score should raise FileNotFoundError if genomes_file is missing."""
    missing_genomes = tmp_path / "missing_genomes.fa"

    with pytest.raises(FileNotFoundError, match="Genomes FASTA file not found"):
        compute_alignment_score(
            genomes_file=missing_genomes,
            atomization_file=mini_geese,
            output_directory=output_dir
        )


# -----------------------------------------------------------------------------
# Test: raises FileNotFoundError if atomization_file is missing
# -----------------------------------------------------------------------------
def test_compute_alignment_score_missing_atomization(mini_fasta: Path, tmp_path: Path, output_dir: Path):
    """compute_alignment_score should raise FileNotFoundError if atomization_file is missing."""
    missing_atomization = tmp_path / "missing_atomization.geese"

    with pytest.raises(FileNotFoundError, match="Atomization file not found"):
        compute_alignment_score(
            genomes_file=mini_fasta,
            atomization_file=missing_atomization,
            output_directory=output_dir
        )


# -----------------------------------------------------------------------------
# Test: raises ValueError for invalid level
# -----------------------------------------------------------------------------
def test_compute_alignment_score_invalid_level(mini_fasta: Path, mini_geese: Path, output_dir: Path, monkeypatch):
    """compute_alignment_score should raise ValueError for an invalid level."""

    monkeypatch.setattr(
        "atomization_scorer.scoring_system.alignment_score.compute_true_alignment",
        lambda *args, **kwargs: mini_geese
    )

    with pytest.raises(ValueError, match="Level must be 'base' or 'interval'"):
        compute_alignment_score(
            genomes_file=mini_fasta,
            atomization_file=mini_geese,
            output_directory=output_dir,
            level="invalid"
        )
