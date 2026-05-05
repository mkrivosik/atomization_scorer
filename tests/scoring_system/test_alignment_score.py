"""
Tests for the alignment_score.py module.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
from pathlib import Path

import pytest

from atomization_scorer.scoring_system import compute_alignment_score


# --------------------------------------------------------------------------------------
# Test: base-level score computation forwards expected arguments
# --------------------------------------------------------------------------------------
def test_compute_alignment_score_base(
    mini_fasta: Path,
    mini_geese: Path,
    output_dir: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    """compute_alignment_score should call base-level metrics with the expected arguments."""
    true_geese = output_dir / "true_atomization" / "true_atomization.geese"
    calls = {}

    def fake_true_alignment(**kwargs):
        calls["true_alignment"] = kwargs
        return true_geese

    def fake_base_metrics(**kwargs):
        calls["base_metrics"] = kwargs
        return 0.85

    def fake_atomization_visualization(**kwargs):
        calls["atomization_visualization"] = kwargs

    monkeypatch.setattr(
        "atomization_scorer.scoring_system.alignment_score.compute_true_alignment",
        fake_true_alignment
    )
    monkeypatch.setattr(
        "atomization_scorer.scoring_system.alignment_score.compute_base_level_metrics",
        fake_base_metrics
    )
    monkeypatch.setattr(
        "atomization_scorer.scoring_system.alignment_score.plot_atomization",
        fake_atomization_visualization
    )

    score = compute_alignment_score(
        genomes_file=mini_fasta,
        atomization_file=mini_geese,
        output_directory=output_dir,
        level="base",
        per_class=False
    )

    assert score == 0.85
    assert calls["true_alignment"] == {
        "genomes_file": mini_fasta,
        "atomization_file": mini_geese,
        "output_directory": output_dir / "true_atomization",
        "mode": "mash",
        "minimum_similarity": 0.95,
        "minimum_alignment_length": 500,
        "minimap2_preset": "asm20",
        "minimap2_secondary_ratio": 0.1,
        "minimap2_emit_cigar": True,
        "run_overlap_diagnostics": False,
        "minimum_report_overlap_length": 0,
        "minimum_plot_overlap_length": 0,
        "overlap_include_reverse": False,
        "run_dotter": True,
        "diagnostics_directory": output_dir / "overlap_diagnostics",
    }
    assert calls["base_metrics"] == {
        "predicted_geese": mini_geese,
        "true_geese": true_geese,
        "output_directory": output_dir / "metrics",
        "per_class": False,
    }
    assert calls["atomization_visualization"] == {
        "genomes_file": mini_fasta,
        "true_atoms_file": true_geese,
        "predicted_atoms_file": mini_geese,
        "output_directory": output_dir / "visualization",
    }


# --------------------------------------------------------------------------------------
# Test: interval-level score computation forwards expected arguments
# --------------------------------------------------------------------------------------
def test_compute_alignment_score_interval(
    mini_fasta: Path,
    mini_geese: Path,
    output_dir: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    """compute_alignment_score should call interval-level metrics with the expected arguments."""
    true_geese = output_dir / "true_atomization" / "true_atomization.geese"
    calls = {}

    def fake_true_alignment(**kwargs):
        calls["true_alignment"] = kwargs
        return true_geese

    def fake_interval_metrics(**kwargs):
        calls["interval_metrics"] = kwargs
        return 0.75

    def fake_atomization_visualization(**kwargs):
        calls["atomization_visualization"] = kwargs

    monkeypatch.setattr(
        "atomization_scorer.scoring_system.alignment_score.compute_true_alignment",
        fake_true_alignment
    )
    monkeypatch.setattr(
        "atomization_scorer.scoring_system.alignment_score.compute_interval_level_metrics",
        fake_interval_metrics
    )
    monkeypatch.setattr(
        "atomization_scorer.scoring_system.alignment_score.plot_atomization",
        fake_atomization_visualization
    )

    score = compute_alignment_score(
        genomes_file=mini_fasta,
        atomization_file=mini_geese,
        output_directory=output_dir,
        level="interval",
        per_class=False,
        minimum_overlap_ratio=0.9
    )

    assert score == 0.75
    assert calls["true_alignment"] == {
        "genomes_file": mini_fasta,
        "atomization_file": mini_geese,
        "output_directory": output_dir / "true_atomization",
        "mode": "mash",
        "minimum_similarity": 0.95,
        "minimum_alignment_length": 500,
        "minimap2_preset": "asm20",
        "minimap2_secondary_ratio": 0.1,
        "minimap2_emit_cigar": True,
        "run_overlap_diagnostics": False,
        "minimum_report_overlap_length": 0,
        "minimum_plot_overlap_length": 0,
        "overlap_include_reverse": False,
        "run_dotter": True,
        "diagnostics_directory": output_dir / "overlap_diagnostics",
    }
    assert calls["interval_metrics"] == {
        "predicted_geese": mini_geese,
        "true_geese": true_geese,
        "output_directory": output_dir / "metrics",
        "per_class": False,
        "minimum_overlap_ratio": 0.9,
    }
    assert calls["atomization_visualization"] == {
        "genomes_file": mini_fasta,
        "true_atoms_file": true_geese,
        "predicted_atoms_file": mini_geese,
        "output_directory": output_dir / "visualization",
    }


# --------------------------------------------------------------------------------------
# Test: non-default pipeline and minimap2 parameters are forwarded
# --------------------------------------------------------------------------------------
def test_compute_alignment_score_forwards_pipeline_parameters(
    mini_fasta: Path,
    mini_geese: Path,
    output_dir: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    """compute_alignment_score should thread non-default pipeline and minimap2 parameters through."""
    true_geese = output_dir / "true_atomization" / "true_atomization.geese"
    calls = {}

    def fake_true_alignment(**kwargs):
        calls["true_alignment"] = kwargs
        return true_geese

    monkeypatch.setattr(
        "atomization_scorer.scoring_system.alignment_score.compute_true_alignment",
        fake_true_alignment
    )
    monkeypatch.setattr(
        "atomization_scorer.scoring_system.alignment_score.compute_interval_level_metrics",
        lambda **_: 0.5
    )
    monkeypatch.setattr(
        "atomization_scorer.scoring_system.alignment_score.plot_atomization",
        lambda **_: None
    )

    compute_alignment_score(
        genomes_file=mini_fasta,
        atomization_file=mini_geese,
        output_directory=output_dir,
        representative_mode="first",
        minimum_similarity=0.8,
        minimum_alignment_length=1000,
        minimap2_preset="asm5",
        minimap2_secondary_ratio=0.05,
        minimap2_emit_cigar=False,
        run_overlap_diagnostics=True,
        minimum_report_overlap_length=100,
        minimum_plot_overlap_length=500,
        overlap_include_reverse=True,
        run_dotter=False,
    )

    assert calls["true_alignment"]["mode"] == "first"
    assert calls["true_alignment"]["minimum_similarity"] == 0.8
    assert calls["true_alignment"]["minimum_alignment_length"] == 1000
    assert calls["true_alignment"]["minimap2_preset"] == "asm5"
    assert calls["true_alignment"]["minimap2_secondary_ratio"] == 0.05
    assert calls["true_alignment"]["minimap2_emit_cigar"] is False
    assert calls["true_alignment"]["run_overlap_diagnostics"] is True
    assert calls["true_alignment"]["minimum_report_overlap_length"] == 100
    assert calls["true_alignment"]["minimum_plot_overlap_length"] == 500
    assert calls["true_alignment"]["overlap_include_reverse"] is True
    assert calls["true_alignment"]["run_dotter"] is False


# --------------------------------------------------------------------------------------
# Test: invalid level does not call gold-standard pipeline
# --------------------------------------------------------------------------------------
def test_compute_alignment_score_invalid_level_does_not_call_true_alignment(
    mini_fasta: Path,
    mini_geese: Path,
    output_dir: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    """compute_alignment_score should reject an invalid level before computing the gold standard."""
    was_called = False

    def fake_true_alignment(**_kwargs):
        nonlocal was_called
        was_called = True
        return mini_geese

    monkeypatch.setattr(
        "atomization_scorer.scoring_system.alignment_score.compute_true_alignment",
        fake_true_alignment
    )

    with pytest.raises(ValueError, match="Level must be 'base' or 'interval'"):
        compute_alignment_score(
            genomes_file=mini_fasta,
            atomization_file=mini_geese,
            output_directory=output_dir,
            level="invalid"
        )

    assert was_called is False


# --------------------------------------------------------------------------------------
# Test: propagates true-alignment failure
# --------------------------------------------------------------------------------------
def test_compute_alignment_score_propagates_true_alignment_failure(
    mini_fasta: Path,
    mini_geese: Path,
    output_dir: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    """compute_alignment_score should propagate failures from gold-standard generation."""
    def fake_true_alignment(**_kwargs):
        raise RuntimeError("true alignment failed")

    monkeypatch.setattr(
        "atomization_scorer.scoring_system.alignment_score.compute_true_alignment",
        fake_true_alignment
    )

    with pytest.raises(RuntimeError, match="true alignment failed"):
        compute_alignment_score(
            genomes_file=mini_fasta,
            atomization_file=mini_geese,
            output_directory=output_dir,
            level="base"
        )


# --------------------------------------------------------------------------------------
# Test: propagates base-level metrics failure
# --------------------------------------------------------------------------------------
def test_compute_alignment_score_propagates_base_metrics_failure(
    mini_fasta: Path,
    mini_geese: Path,
    output_dir: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    """compute_alignment_score should propagate failures from base-level metrics."""
    monkeypatch.setattr(
        "atomization_scorer.scoring_system.alignment_score.compute_true_alignment",
        lambda **_kwargs: mini_geese
    )

    def fake_base_metrics(**_kwargs):
        raise RuntimeError("base metrics failed")

    monkeypatch.setattr(
        "atomization_scorer.scoring_system.alignment_score.compute_base_level_metrics",
        fake_base_metrics
    )
    monkeypatch.setattr(
        "atomization_scorer.scoring_system.alignment_score.plot_atomization",
        lambda **_kwargs: None
    )

    with pytest.raises(RuntimeError, match="base metrics failed"):
        compute_alignment_score(
            genomes_file=mini_fasta,
            atomization_file=mini_geese,
            output_directory=output_dir,
            level="base"
        )


# --------------------------------------------------------------------------------------
# Test: propagates interval-level metrics failure
# --------------------------------------------------------------------------------------
def test_compute_alignment_score_propagates_interval_metrics_failure(
    mini_fasta: Path,
    mini_geese: Path,
    output_dir: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    """compute_alignment_score should propagate failures from interval-level metrics."""
    monkeypatch.setattr(
        "atomization_scorer.scoring_system.alignment_score.compute_true_alignment",
        lambda **_kwargs: mini_geese
    )

    def fake_interval_metrics(**_kwargs):
        raise RuntimeError("interval metrics failed")

    monkeypatch.setattr(
        "atomization_scorer.scoring_system.alignment_score.compute_interval_level_metrics",
        fake_interval_metrics
    )
    monkeypatch.setattr(
        "atomization_scorer.scoring_system.alignment_score.plot_atomization",
        lambda **_kwargs: None
    )

    with pytest.raises(RuntimeError, match="interval metrics failed"):
        compute_alignment_score(
            genomes_file=mini_fasta,
            atomization_file=mini_geese,
            output_directory=output_dir,
            level="interval"
        )


# --------------------------------------------------------------------------------------
# Test: raises FileNotFoundError if an input file is missing
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("missing_side", "expected_message"),
    [
        ("genomes", "Genomes FASTA file not found"),
        ("atomization", "Atomization file not found"),
    ],
)
def test_compute_alignment_score_missing_input_file(
    mini_fasta: Path,
    mini_geese: Path,
    tmp_path: Path,
    output_dir: Path,
    missing_side: str,
    expected_message: str,
):
    """compute_alignment_score should raise FileNotFoundError if an input file is missing."""
    genomes_file = mini_fasta
    atomization_file = mini_geese

    if missing_side == "genomes":
        genomes_file = tmp_path / "missing_genomes.fa"
    else:
        atomization_file = tmp_path / "missing_atomization.geese"

    with pytest.raises(FileNotFoundError, match=expected_message):
        compute_alignment_score(
            genomes_file=genomes_file,
            atomization_file=atomization_file,
            output_directory=output_dir
        )
