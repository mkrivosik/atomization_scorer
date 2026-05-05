"""
Tests for the overall_score.py module.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
from pathlib import Path
from unittest.mock import Mock

import pytest

from atomization_scorer.scoring_system import compute_overall_score


# --------------------------------------------------------------------------------------
# Test: computes overall score and forwards alignment parameters
# --------------------------------------------------------------------------------------
def test_compute_overall_score(
    mini_fasta: Path,
    mini_geese: Path,
    output_dir: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    """compute_overall_score should forward arguments and combine alignment and coverage scores."""
    mock_compute_alignment_score = Mock(return_value=0.8)
    mock_compute_coverage_score = Mock(return_value=0.9)

    monkeypatch.setattr(
        "atomization_scorer.scoring_system.overall_score.compute_alignment_score",
        mock_compute_alignment_score
    )
    monkeypatch.setattr(
        "atomization_scorer.scoring_system.overall_score.compute_coverage_score",
        mock_compute_coverage_score
    )

    overall_score = compute_overall_score(
        genomes_file=mini_fasta,
        atomization_file=mini_geese,
        output_directory=output_dir,
        level="base",
        per_class=False,
        minimum_overlap_ratio=0.65,
        alignment_weight=0.6,
        coverage_weight=0.4
    )

    mock_compute_alignment_score.assert_called_once_with(
        genomes_file=mini_fasta,
        atomization_file=mini_geese,
        output_directory=output_dir,
        level="base",
        per_class=False,
        minimum_overlap_ratio=0.65,
        representative_mode="mash",
        minimum_similarity=0.95,
        minimum_alignment_length=500,
        minimap2_preset="asm20",
        minimap2_secondary_ratio=0.1,
        minimap2_emit_cigar=True,
        run_overlap_diagnostics=False,
        minimum_report_overlap_length=0,
        minimum_plot_overlap_length=0,
        overlap_include_reverse=False,
        run_dotter=True,
    )
    mock_compute_coverage_score.assert_called_once_with(
        genomes_file=mini_fasta,
        atomization_file=mini_geese,
        per_class=False,
    )
    assert overall_score == (0.8 ** 0.6) * (0.9 ** 0.4)


# --------------------------------------------------------------------------------------
# Test: per-class overall score combines per-class alignment and coverage
# --------------------------------------------------------------------------------------
def test_compute_overall_score_per_class(
    mini_fasta: Path,
    mini_geese: Path,
    output_dir: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    """compute_overall_score with per_class=True should return per-class weighted geometric means."""
    alignment_results = [
        {"Class": 1, "F1-score": 0.8},
        {"Class": 2, "F1-score": 0.6},
    ]
    coverage_results = [
        {"Class": 1, "Coverage": 0.5},
        {"Class": 2, "Coverage": 0.3},
    ]

    monkeypatch.setattr(
        "atomization_scorer.scoring_system.overall_score.compute_alignment_score",
        lambda **_: alignment_results,
    )
    monkeypatch.setattr(
        "atomization_scorer.scoring_system.overall_score.compute_coverage_score",
        lambda **_: coverage_results,
    )

    result = compute_overall_score(
        genomes_file=mini_fasta,
        atomization_file=mini_geese,
        output_directory=output_dir,
        per_class=True,
        alignment_weight=0.7,
        coverage_weight=0.3,
    )

    assert isinstance(result, list)
    assert len(result) == 2

    class1 = next(r for r in result if r["Class"] == 1)
    class2 = next(r for r in result if r["Class"] == 2)

    assert class1["Alignment_score"] == 0.8
    assert class1["Coverage_score"] == 0.5
    assert class1["Overall"] == pytest.approx((0.8 ** 0.7) * (0.5 ** 0.3))

    assert class2["Alignment_score"] == 0.6
    assert class2["Coverage_score"] == 0.3
    assert class2["Overall"] == pytest.approx((0.6 ** 0.7) * (0.3 ** 0.3))


# --------------------------------------------------------------------------------------
# Test: per-class fills missing classes with 0.0
# --------------------------------------------------------------------------------------
def test_compute_overall_score_per_class_missing_class(
    mini_fasta: Path,
    mini_geese: Path,
    output_dir: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    """compute_overall_score per-class should default to 0.0 for classes missing from one side."""
    monkeypatch.setattr(
        "atomization_scorer.scoring_system.overall_score.compute_alignment_score",
        lambda **_: [{"Class": 1, "F1-score": 0.8}, {"Class": 3, "F1-score": 0.5}],
    )
    monkeypatch.setattr(
        "atomization_scorer.scoring_system.overall_score.compute_coverage_score",
        lambda **_: [{"Class": 1, "Coverage": 0.4}, {"Class": 2, "Coverage": 0.2}],
    )

    result = compute_overall_score(
        genomes_file=mini_fasta,
        atomization_file=mini_geese,
        output_directory=output_dir,
        per_class=True,
    )

    assert isinstance(result, list)
    classes = {r["Class"] for r in result}
    assert classes == {1, 2, 3}

    class2 = next(r for r in result if r["Class"] == 2)
    assert class2["Alignment_score"] == 0.0

    class3 = next(r for r in result if r["Class"] == 3)
    assert class3["Coverage_score"] == 0.0


# --------------------------------------------------------------------------------------
# Test: returns expected edge-case score
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("alignment_score", "coverage_score", "expected_score"),
    [
        (0.5, 0.0, 0.0),
        (0.0, 0.5, 0.0),
        (1.0, 1.0, 1.0),
    ]
)
def test_compute_overall_score_edge_cases(
    mini_fasta: Path,
    mini_geese: Path,
    output_dir: Path,
    monkeypatch: pytest.MonkeyPatch,
    alignment_score: float,
    coverage_score: float,
    expected_score: float,
):
    """compute_overall_score should return the expected score for edge-case inputs."""
    monkeypatch.setattr(
        "atomization_scorer.scoring_system.overall_score.compute_alignment_score",
        lambda **_kwargs: alignment_score
    )
    monkeypatch.setattr(
        "atomization_scorer.scoring_system.overall_score.compute_coverage_score",
        lambda **_kwargs: coverage_score
    )

    overall = compute_overall_score(
        genomes_file=mini_fasta,
        atomization_file=mini_geese,
        output_directory=output_dir
    )
    assert overall == expected_score


# --------------------------------------------------------------------------------------
# Test: output_directory is created if it does not exist
# --------------------------------------------------------------------------------------
def test_compute_overall_score_creates_output_dir(
    mini_fasta: Path,
    mini_geese: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
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


# --------------------------------------------------------------------------------------
# Test: raises ValueError for invalid weights
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("alignment_weight", "coverage_weight", "expected_message"),
    [
        (-0.1, 1.1, "Alignment and coverage weights must be non-negative."),
        (0.7, -0.3, "Alignment and coverage weights must be non-negative."),
        (0.6, 0.3, "Alignment and coverage weights must sum to 1.0."),
    ]
)
def test_compute_overall_score_invalid_weights(
    mini_fasta: Path,
    mini_geese: Path,
    output_dir: Path,
    monkeypatch: pytest.MonkeyPatch,
    alignment_weight: float,
    coverage_weight: float,
    expected_message: str,
):
    """compute_overall_score should reject invalid geometric-mean weights."""
    mock_compute_alignment_score = Mock()
    mock_compute_coverage_score = Mock()

    monkeypatch.setattr(
        "atomization_scorer.scoring_system.overall_score.compute_alignment_score",
        mock_compute_alignment_score
    )
    monkeypatch.setattr(
        "atomization_scorer.scoring_system.overall_score.compute_coverage_score",
        mock_compute_coverage_score
    )

    with pytest.raises(ValueError, match=expected_message):
        compute_overall_score(
            genomes_file=mini_fasta,
            atomization_file=mini_geese,
            output_directory=output_dir,
            alignment_weight=alignment_weight,
            coverage_weight=coverage_weight
        )

    mock_compute_alignment_score.assert_not_called()
    mock_compute_coverage_score.assert_not_called()


# --------------------------------------------------------------------------------------
# Test: propagates alignment-score failure
# --------------------------------------------------------------------------------------
def test_compute_overall_score_propagates_alignment_score_failure(
    mini_fasta: Path,
    mini_geese: Path,
    output_dir: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    """compute_overall_score should propagate compute_alignment_score failures."""
    mock_compute_coverage_score = Mock()

    def fake_compute_alignment_score(**_kwargs):
        raise RuntimeError("alignment score failed")

    monkeypatch.setattr(
        "atomization_scorer.scoring_system.overall_score.compute_alignment_score",
        fake_compute_alignment_score
    )
    monkeypatch.setattr(
        "atomization_scorer.scoring_system.overall_score.compute_coverage_score",
        mock_compute_coverage_score
    )

    with pytest.raises(RuntimeError, match="alignment score failed"):
        compute_overall_score(
            genomes_file=mini_fasta,
            atomization_file=mini_geese,
            output_directory=output_dir
        )

    mock_compute_coverage_score.assert_not_called()


# --------------------------------------------------------------------------------------
# Test: propagates coverage-score failure
# --------------------------------------------------------------------------------------
def test_compute_overall_score_propagates_coverage_score_failure(
    mini_fasta: Path,
    mini_geese: Path,
    output_dir: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    """compute_overall_score should propagate compute_coverage_score failures."""
    monkeypatch.setattr(
        "atomization_scorer.scoring_system.overall_score.compute_alignment_score",
        lambda **_kwargs: 0.8
    )

    def fake_compute_coverage_score(**_kwargs):
        raise RuntimeError("coverage score failed")

    monkeypatch.setattr(
        "atomization_scorer.scoring_system.overall_score.compute_coverage_score",
        fake_compute_coverage_score
    )

    with pytest.raises(RuntimeError, match="coverage score failed"):
        compute_overall_score(
            genomes_file=mini_fasta,
            atomization_file=mini_geese,
            output_directory=output_dir
        )


# --------------------------------------------------------------------------------------
# Test: raises FileNotFoundError if an input file is missing
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("missing_input", "expected_message"),
    [
        ("genomes", "Genomes FASTA file not found"),
        ("atomization", "Atomization file not found"),
    ]
)
def test_compute_overall_score_missing_input_file(
    mini_fasta: Path,
    mini_geese: Path,
    output_dir: Path,
    missing_input: str,
    expected_message: str
):
    """compute_overall_score should raise FileNotFoundError if an input file does not exist."""
    genomes_file = output_dir / "not_exist.fa" if missing_input == "genomes" else mini_fasta
    atomization_file = output_dir / "not_exist.geese" if missing_input == "atomization" else mini_geese

    with pytest.raises(FileNotFoundError, match=expected_message):
        compute_overall_score(
            genomes_file=genomes_file,
            atomization_file=atomization_file,
            output_directory=output_dir
        )
