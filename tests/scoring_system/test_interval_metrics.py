"""
Tests for the interval_metrics.py module.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
from pathlib import Path

import pandas as pd

from atomization_scorer.scoring_system import compute_interval_level_metrics


# --------------------------------------------------------------------------------------
# Helper: create minimal predicted and true GEESE files
# --------------------------------------------------------------------------------------
def create_minimal_geese(tmp_path: Path):
    """Create minimal predicted and true GEESE files for interval-level testing."""
    predicted_geese = tmp_path / "predicted.geese"
    predicted_geese.write_text(
        "#name\tatom_nr\tclass\tstrand\tstart\tend\n"
        "sequence1\t1\t1\t+\t0\t10\n"   # class = 1: TP = 1, FP = 2, FN = 1
        "sequence1\t2\t1\t+\t20\t30\n"  
        "sequence2\t3\t1\t+\t0\t10\n"
    )

    true_geese = tmp_path / "true.geese"
    true_geese.write_text(
        "#name\tatom_nr\tclass\tstrand\tstart\tend\n"
        "sequence1\t1\t1\t+\t1\t11\n"
        "sequence2\t2\t2\t+\t0\t20\n"   # class = 2: TP = 0, FP = 0, FN = 1
        "sequence2\t3\t1\t+\t8\t20\n"
        ""
    )

    return predicted_geese, true_geese


# --------------------------------------------------------------------------------------
# Helper: create GEESE files for class-mismatch testing
# --------------------------------------------------------------------------------------
def create_class_mismatch_geese(tmp_path: Path):
    """Create predicted and true GEESE files where overlap exists but classes differ."""
    predicted_geese = tmp_path / "predicted_class_mismatch.geese"
    predicted_geese.write_text(
        "#name\tatom_nr\tclass\tstrand\tstart\tend\n"
        "sequence1\t1\t1\t+\t0\t10\n"
    )

    true_geese = tmp_path / "true_class_mismatch.geese"
    true_geese.write_text(
        "#name\tatom_nr\tclass\tstrand\tstart\tend\n"
        "sequence1\t1\t2\t+\t0\t10\n"
    )

    return predicted_geese, true_geese


# --------------------------------------------------------------------------------------
# Helper: create GEESE files where one true interval can only be matched once
# --------------------------------------------------------------------------------------
def create_single_true_match_geese(tmp_path: Path):
    """Create predicted and true GEESE files where two predictions compete for one true interval."""
    predicted_geese = tmp_path / "predicted_single_true_match.geese"
    predicted_geese.write_text(
        "#name\tatom_nr\tclass\tstrand\tstart\tend\n"
        "sequence1\t1\t1\t+\t0\t10\n"
        "sequence1\t2\t1\t+\t0\t10\n"
    )

    true_geese = tmp_path / "true_single_true_match.geese"
    true_geese.write_text(
        "#name\tatom_nr\tclass\tstrand\tstart\tend\n"
        "sequence1\t1\t1\t+\t0\t10\n"
    )

    return predicted_geese, true_geese


# --------------------------------------------------------------------------------------
# Helper: create GEESE files for exact min_overlap_ratio testing
# --------------------------------------------------------------------------------------
def create_min_overlap_ratio_geese(tmp_path: Path):
    """Create predicted and true GEESE files where one match depends on the overlap threshold."""
    predicted_geese = tmp_path / "predicted_min_overlap.geese"
    predicted_geese.write_text(
        "#name\tatom_nr\tclass\tstrand\tstart\tend\n"
        "sequence1\t1\t1\t+\t0\t10\n"
    )

    true_geese = tmp_path / "true_min_overlap.geese"
    true_geese.write_text(
        "#name\tatom_nr\tclass\tstrand\tstart\tend\n"
        "sequence1\t1\t1\t+\t2\t10\n"
    )

    return predicted_geese, true_geese


# --------------------------------------------------------------------------------------
# Test: basic interval-level scoring (overall)
# --------------------------------------------------------------------------------------
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

    predicted_file = output_dir / "interval_predicted_statuses.tsv"
    true_file = output_dir / "interval_true_statuses.tsv"
    out_file = output_dir / "interval_metrics_overall.tsv"

    assert predicted_file.exists()
    assert true_file.exists()
    assert out_file.is_file()


# --------------------------------------------------------------------------------------
# Test: per-class interval-level scoring
# --------------------------------------------------------------------------------------
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

    predicted_file = output_dir / "interval_predicted_statuses.tsv"
    true_file = output_dir / "interval_true_statuses.tsv"
    out_file = output_dir / "interval_metrics_per_class.tsv"

    assert predicted_file.exists()
    assert true_file.exists()
    assert out_file.is_file()


# --------------------------------------------------------------------------------------
# Test: interval status files are written
# --------------------------------------------------------------------------------------
def test_interval_level_status_files_created(tmp_path: Path, output_dir: Path):
    """compute_interval_level_metrics should create predicted and true status TSV files."""
    predicted_geese, true_geese = create_minimal_geese(tmp_path)

    compute_interval_level_metrics(
        predicted_geese=predicted_geese,
        true_geese=true_geese,
        output_directory=output_dir
    )

    predicted_status = output_dir / "interval_predicted_statuses.tsv"
    true_status = output_dir / "interval_true_statuses.tsv"

    assert predicted_status.is_file()
    assert true_status.is_file()

    pred_df = pd.read_csv(predicted_status, sep="\t")
    true_df = pd.read_csv(true_status, sep="\t")

    assert "status" in pred_df.columns
    assert "status" in true_df.columns
    assert list(zip(pred_df["name"], pred_df["start"], pred_df["end"], pred_df["status"])) == [
        ("sequence1", 0, 10, "TP"),
        ("sequence2", 0, 10, "FP"),
        ("sequence1", 20, 30, "FP"),
    ]
    assert list(zip(true_df["name"], true_df["start"], true_df["end"], true_df["status"])) == [
        ("sequence2", 0, 20, "FN"),
        ("sequence1", 1, 11, "TP"),
        ("sequence2", 8, 20, "FN"),
    ]


# --------------------------------------------------------------------------------------
# Test: min_overlap_ratio filters matches
# --------------------------------------------------------------------------------------
def test_interval_level_min_overlap_ratio_effect(tmp_path: Path, output_dir: Path):
    """compute_interval_level_metrics should change the outcome when overlap falls below the threshold."""
    predicted_geese, true_geese = create_min_overlap_ratio_geese(tmp_path)

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
        min_overlap_ratio=0.9
    )

    assert score_low == 1.0
    assert score_high == 0.0


# --------------------------------------------------------------------------------------
# Test: class mismatch does not count as a match
# --------------------------------------------------------------------------------------
def test_interval_level_class_mismatch(tmp_path: Path, output_dir: Path):
    """compute_interval_level_metrics should reject overlaps when predicted and true classes differ."""
    predicted_geese, true_geese = create_class_mismatch_geese(tmp_path)

    score = compute_interval_level_metrics(
        predicted_geese=predicted_geese,
        true_geese=true_geese,
        output_directory=output_dir,
        min_overlap_ratio=0.8
    )

    predicted_status = pd.read_csv(output_dir / "interval_predicted_statuses.tsv", sep="\t")
    true_status = pd.read_csv(output_dir / "interval_true_statuses.tsv", sep="\t")

    assert score == 0.0
    assert predicted_status["status"].tolist() == ["FP"]
    assert true_status["status"].tolist() == ["FN"]


# --------------------------------------------------------------------------------------
# Test: one true interval can be matched only once
# --------------------------------------------------------------------------------------
def test_interval_level_true_interval_matched_only_once(tmp_path: Path, output_dir: Path):
    """compute_interval_level_metrics should allow only one predicted interval to match a true interval."""
    predicted_geese, true_geese = create_single_true_match_geese(tmp_path)

    score = compute_interval_level_metrics(
        predicted_geese=predicted_geese,
        true_geese=true_geese,
        output_directory=output_dir,
        min_overlap_ratio=0.8
    )

    predicted_status = pd.read_csv(output_dir / "interval_predicted_statuses.tsv", sep="\t")
    true_status = pd.read_csv(output_dir / "interval_true_statuses.tsv", sep="\t")

    assert score == 2 / 3
    assert predicted_status["status"].tolist() == ["TP", "FP"]
    assert true_status["status"].tolist() == ["TP"]
