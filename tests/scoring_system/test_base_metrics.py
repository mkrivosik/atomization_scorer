"""
Tests for the interval_metrics.py functions.
"""

import pandas as pd
from pathlib import Path
from atomization_scorer.scoring_system import compute_base_level_metrics


# ---------------------------------------------------------------------
# Helper: create minimal predicted and true GEESE files
# ---------------------------------------------------------------------
def create_minimal_geese(tmp_path: Path):
    """Create minimal predicted and true GEESE files for base-level testing."""
    predicted_geese = tmp_path / "predicted.geese"
    predicted_geese.write_text(
        "#name\tatom_nr\tclass\tstrand\tstart\tend\n"
        "sequence1\t1\t1\t+\t0\t19\n"   # class = 1: TP = 10, FP = 10, FN = 30
        "sequence2\t2\t2\t+\t10\t29\n"  # class = 2: TP = 0, FP = 20, FN = 0
    )

    true_geese = tmp_path / "true.geese"
    true_geese.write_text(
        "#name\tatom_nr\tclass\tstrand\tstart\tend\n"
        "sequence1\t1\t1\t+\t10\t29\n"
        "sequence2\t2\t1\t+\t0\t19\n"
    )

    return predicted_geese, true_geese


# ---------------------------------------------------------------------
# Test: basic base-level scoring (overall)
# ---------------------------------------------------------------------
def test_compute_base_level_metrics(tmp_path: Path, output_dir: Path):
    """compute_base_level_metrics should compute overall base-level metrics."""
    predicted_geese, true_geese = create_minimal_geese(tmp_path)

    metrics = compute_base_level_metrics(
        predicted_geese=predicted_geese,
        true_geese=true_geese,
        output_directory=output_dir,
        per_class=False
    )

    assert metrics == 0.25

    predicted_file = output_dir / "base_predicted_status.tsv"
    true_file = output_dir / "base_true_status.tsv"
    out_file = output_dir / "base_metrics_overall.tsv"

    assert predicted_file.exists()
    assert true_file.exists()
    assert out_file.is_file()


# ---------------------------------------------------------------------
# Test: basic base-level scoring (overall)
# ---------------------------------------------------------------------
def test_compute_per_class_level_metrics(tmp_path: Path, output_dir: Path):
    """compute_base_level_metrics should compute base-level metrics per class."""
    predicted_geese, true_geese = create_minimal_geese(tmp_path)

    metrics = compute_base_level_metrics(
        predicted_geese=predicted_geese,
        true_geese=true_geese,
        output_directory=output_dir,
        per_class=True
    )

    assert isinstance(metrics, list)
    assert len(metrics) == 2
    assert all("Class" in m and "F1-score" in m for m in metrics)
    assert metrics[0]["Class"] == 1
    assert metrics[0]["F1-score"] == 1/3  # Based on TP = 10, FP = 10, FN = 30
    assert metrics[1]["Class"] == 2
    assert metrics[1]["F1-score"] == 0.0  # Based on TP = 0, FP = 20, FN = 0

    predicted_file = output_dir / "base_predicted_status.tsv"
    true_file = output_dir / "base_true_status.tsv"
    out_file = output_dir / "base_metrics_per_class.tsv"

    assert predicted_file.exists()
    assert true_file.exists()
    assert out_file.is_file()


# ---------------------------------------------------------------------
# Test: empty GEESE files
# ---------------------------------------------------------------------
def test_empty_geese(tmp_path: Path, output_dir: Path):
    """compute_base_level_metrics should correctly calculate 0.0 TP, FP, FN."""
    empty_df = pd.DataFrame(columns=["name", "atom_nr", "class","strand", "start", "end"])

    predicted_geese = tmp_path / "predicted_empty.geese"
    true_geese = tmp_path / "true_empty.geese"

    empty_df.to_csv(predicted_geese, sep="\t", index=False)
    empty_df.to_csv(true_geese, sep="\t", index=False)

    metrics = compute_base_level_metrics(
        predicted_geese=predicted_geese,
        true_geese=true_geese,
        output_directory=output_dir,
        per_class=False
    )

    assert metrics == 0.0


# ---------------------------------------------------------------------
# Test: partial class overlap with same classes
# ---------------------------------------------------------------------
def test_partial_class_overlap_same_classes(tmp_path: Path, output_dir: Path):
    """
    compute_base_level_metrics should compute metrics when prediction and true values overlap within the same class.
    """
    predicted_geese = tmp_path / "predicted_partial_overlap.geese"
    predicted_geese.write_text(
        "#name\tatom_nr\tclass\tstrand\tstart\tend\n"
        "sequence1\t1\t1\t+\t0\t19\n"   # TP = 10, FP = 10, FN = 10
    )

    true_geese = tmp_path / "true_partial_overlap.geese"
    true_geese.write_text(
        "#name\tatom_nr\tclass\tstrand\tstart\tend\n"
        "sequence1\t1\t1\t+\t10\t29\n"
    )

    metrics = compute_base_level_metrics(
        predicted_geese=predicted_geese,
        true_geese=true_geese,
        output_directory=output_dir,
        per_class=False
    )

    assert metrics == 0.5


# ---------------------------------------------------------------------
# Test: partial class overlap
# ---------------------------------------------------------------------
def test_partial_class_overlap_different_classes(tmp_path: Path, output_dir: Path):
    """
    compute_base_level_metrics should compute metrics when prediction and true values overlap across different classes.
    """
    predicted_geese = tmp_path / "predicted_partial_overlap.geese"
    predicted_geese.write_text(
        "#name\tatom_nr\tclass\tstrand\tstart\tend\n"
        "sequence1\t1\t1\t+\t0\t19\n"   # TP = 0, FP = 20, FN = 20
    )

    true_geese = tmp_path / "true_partial_overlap.geese"
    true_geese.write_text(
        "#name\tatom_nr\tclass\tstrand\tstart\tend\n"
        "sequence1\t1\t2\t+\t10\t29\n"
    )

    metrics = compute_base_level_metrics(
        predicted_geese=predicted_geese,
        true_geese=true_geese,
        output_directory=output_dir,
        per_class=False
    )

    assert metrics == 0.0


# ---------------------------------------------------------------------
# Test: full class overlap
# ---------------------------------------------------------------------
def test_full_class_overlap(tmp_path: Path, output_dir: Path):
    """compute_base_level_metrics should return 1.0 when prediction and true values are identical."""
    predicted_geese = tmp_path / "predicted_full_overlap.geese"
    predicted_geese.write_text(
        "#name\tatom_nr\tclass\tstrand\tstart\tend\n"
        "sequence1\t1\t1\t+\t0\t9\n"
        "sequence1\t2\t1\t+\t10\t19\n"
    )

    true_geese = tmp_path / "true_full_overlap.geese"
    true_geese.write_text(
        "#name\tatom_nr\tclass\tstrand\tstart\tend\n"
        "sequence1\t1\t1\t+\t0\t9\n"
        "sequence1\t2\t1\t+\t10\t19\n"
    )

    metrics = compute_base_level_metrics(
        predicted_geese=predicted_geese,
        true_geese=true_geese,
        output_directory=output_dir,
        per_class=False
    )

    assert metrics == 1.0
