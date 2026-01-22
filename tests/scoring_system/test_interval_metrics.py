"""
Tests for the interval_metrics.py functions.
"""

from pathlib import Path

import pandas as pd

from atomization_scorer.scoring_system import compute_interval_level_metrics


# ---------------------------------------------------------------------
# Helper: create minimal predicted and true GEESE files
# ---------------------------------------------------------------------
def create_minimal_geese(tmp_path: Path):
    """Create minimal predicted and true GEESE files for interval-level testing."""
    predicted_geese = tmp_path / "predicted.geese"
    predicted_geese.write_text(
        "#name\tatom_nr\tclass\tstrand\tstart\tend\n"
        "sequence1\t1\t1\t+\t0\t9\n"    # class = 1: TP = 1, FP = 2, FN = 1
        "sequence1\t2\t1\t+\t20\t30\n"  
        "sequence2\t3\t1\t+\t0\t9\n"
    )

    true_geese = tmp_path / "true.geese"
    true_geese.write_text(
        "#name\tatom_nr\tclass\tstrand\tstart\tend\n"
        "sequence1\t1\t1\t+\t1\t10\n"
        "sequence2\t2\t2\t+\t0\t20\n"   # class = 2: TP = 0, FP = 0, FN = 1
        "sequence2\t3\t1\t+\t8\t20\n"
        ""
    )

    return predicted_geese, true_geese


# ---------------------------------------------------------------------
# Test: basic interval-level scoring (overall)
# ---------------------------------------------------------------------
def test_compute_interval_level_metrics_overall(tmp_path: Path, output_dir: Path):
    """compute_interval_level_metrics should compute overall interval-level metrics."""
    predicted_geese, true_geese = create_minimal_geese(tmp_path)

    score = compute_interval_level_metrics(
        predicted_geese=predicted_geese,
        true_geese=true_geese,
        output_directory=output_dir,
        per_class=False,
        min_overlap_ratio=0.8
    )

    assert score == 1/3

    predicted_file = output_dir / "interval_predicted_status.tsv"
    true_file = output_dir / "interval_true_status.tsv"
    out_file = output_dir / "interval_metrics_overall.tsv"

    assert predicted_file.exists()
    assert true_file.exists()
    assert out_file.is_file()


# ---------------------------------------------------------------------
# Test: per-class interval-level scoring
# ---------------------------------------------------------------------
def test_compute_interval_level_metrics_per_class(tmp_path: Path, output_dir: Path):
    """compute_interval_level_metrics should compute interval-level metrics per class."""
    predicted_geese, true_geese = create_minimal_geese(tmp_path)

    metrics = compute_interval_level_metrics(
        predicted_geese=predicted_geese,
        true_geese=true_geese,
        output_directory=output_dir,
        per_class=True,
        min_overlap_ratio=0.8
    )

    assert isinstance(metrics, list)
    assert len(metrics) == 2
    assert all("Class" in m and "F1-score" in m for m in metrics)
    assert metrics[0]["Class"] == 1
    assert metrics[0]["F1-score"] == 0.4    # Based on TP = 1, FP = 2, FN = 1
    assert metrics[1]["Class"] == 2
    assert metrics[1]["F1-score"] == 0.0    # Based on TP = 0, FP = 0, FN = 1

    predicted_file = output_dir / "interval_predicted_status.tsv"
    true_file = output_dir / "interval_true_status.tsv"
    out_file = output_dir / "interval_metrics_per_class.tsv"

    assert predicted_file.exists()
    assert true_file.exists()
    assert out_file.is_file()


# ---------------------------------------------------------------------
# Test: interval status files are written
# ---------------------------------------------------------------------
def test_interval_level_status_files_created(tmp_path: Path, output_dir: Path):
    """compute_interval_level_metrics should create predicted and true status TSV files."""
    predicted_geese, true_geese = create_minimal_geese(tmp_path)

    compute_interval_level_metrics(
        predicted_geese=predicted_geese,
        true_geese=true_geese,
        output_directory=output_dir
    )

    predicted_status = output_dir / "interval_predicted_status.tsv"
    true_status = output_dir / "interval_true_status.tsv"

    assert predicted_status.is_file()
    assert true_status.is_file()

    pred_df = pd.read_csv(predicted_status, sep="\t")
    true_df = pd.read_csv(true_status, sep="\t")

    assert "status" in pred_df.columns
    assert "status" in true_df.columns
    assert set(pred_df["status"]).issubset({"TP", "FP"})
    assert set(true_df["status"]).issubset({"TP", "FN"})


# ---------------------------------------------------------------------
# Test: min_overlap_ratio filters matches
# ---------------------------------------------------------------------
def test_interval_level_min_overlap_ratio_effect(tmp_path: Path, output_dir: Path):
    """compute_interval_level_metrics should reduce True Positives with higher min_overlap_ratio."""
    predicted_geese, true_geese = create_minimal_geese(tmp_path)

    score_low = compute_interval_level_metrics(
        predicted_geese=predicted_geese,
        true_geese=true_geese,
        output_directory=output_dir,
        min_overlap_ratio=0.8
    )

    score_high = compute_interval_level_metrics(
        predicted_geese=predicted_geese,
        true_geese=true_geese,
        output_directory=output_dir,
        min_overlap_ratio=1.0
    )

    assert score_high <= score_low
