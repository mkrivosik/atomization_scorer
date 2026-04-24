"""
Tests for the base_metrics.py functions.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
from pathlib import Path

import pandas as pd
import pytest

from atomization_scorer.scoring_system import compute_base_level_metrics


# --------------------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------------------
def create_minimal_geese(tmp_path: Path) -> tuple[Path, Path]:
    """Create minimal predicted and true GEESE files for base-level testing."""
    predicted_geese = tmp_path / "predicted.geese"
    predicted_geese.write_text(
        "#name\tatom_nr\tclass\tstrand\tstart\tend\n"
        "sequence1\t1\t1\t+\t0\t20\n"
        "sequence2\t2\t2\t+\t10\t30\n"
    )

    true_geese = tmp_path / "true.geese"
    true_geese.write_text(
        "#name\tatom_nr\tclass\tstrand\tstart\tend\n"
        "sequence1\t1\t1\t+\t10\t30\n"
        "sequence2\t2\t1\t+\t0\t20\n"
    )

    return predicted_geese, true_geese


def create_partial_overlap_geese(tmp_path: Path) -> tuple[Path, Path]:
    """Create a simple same-class partial overlap case with TP, FP, and FN fragments."""
    predicted_geese = tmp_path / "predicted_partial_overlap.geese"
    predicted_geese.write_text(
        "#name\tatom_nr\tclass\tstrand\tstart\tend\n"
        "sequence1\t1\t1\t+\t0\t20\n"
    )

    true_geese = tmp_path / "true_partial_overlap.geese"
    true_geese.write_text(
        "#name\tatom_nr\tclass\tstrand\tstart\tend\n"
        "sequence1\t1\t1\t+\t10\t30\n"
    )

    return predicted_geese, true_geese


# --------------------------------------------------------------------------------------
# Test: raises FileNotFoundError if an input file is missing
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("missing_side", "expected_message"),
    [
        ("predicted", "Predicted GEESE file not found"),
        ("true", "True GEESE file not found"),
    ],
)
def test_compute_base_level_metrics_missing_input_file(
    tmp_path: Path,
    output_dir: Path,
    missing_side: str,
    expected_message: str,
):
    """compute_base_level_metrics should raise FileNotFoundError if an input file is missing."""
    predicted_geese, true_geese = create_minimal_geese(tmp_path)

    if missing_side == "predicted":
        predicted_geese = tmp_path / "missing_predicted.geese"
    else:
        true_geese = tmp_path / "missing_true.geese"

    with pytest.raises(FileNotFoundError, match=expected_message):
        compute_base_level_metrics(
            predicted_geese=predicted_geese,
            true_geese=true_geese,
            output_directory=output_dir,
            per_class=False
        )


# --------------------------------------------------------------------------------------
# Test: raises ValueError for invalid intervals
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("filename", "row", "expected_message"),
    [
        (
            "predicted_zero_length.geese",
            "sequence1\t1\t1\t+\t10\t10\n",
            r"Invalid predicted interval for sequence 'sequence1', class 1: start=10, end=10",
        ),
        (
            "predicted_reversed.geese",
            "sequence1\t1\t1\t+\t20\t10\n",
            r"Invalid predicted interval for sequence 'sequence1', class 1: start=20, end=10",
        ),
    ],
)
def test_compute_base_level_metrics_invalid_predicted_interval(
    tmp_path: Path,
    output_dir: Path,
    filename: str,
    row: str,
    expected_message: str,
):
    """compute_base_level_metrics should raise ValueError for invalid predicted intervals."""
    predicted_geese = tmp_path / filename
    predicted_geese.write_text("#name\tatom_nr\tclass\tstrand\tstart\tend\n" + row)

    true_geese = tmp_path / "valid_true.geese"
    true_geese.write_text(
        "#name\tatom_nr\tclass\tstrand\tstart\tend\n"
        "sequence1\t1\t1\t+\t0\t20\n"
    )

    with pytest.raises(ValueError, match=expected_message):
        compute_base_level_metrics(
            predicted_geese=predicted_geese,
            true_geese=true_geese,
            output_directory=output_dir,
            per_class=False
        )


# --------------------------------------------------------------------------------------
# Test: touching intervals do not overlap at the shared boundary
# --------------------------------------------------------------------------------------
def test_touching_intervals_do_not_overlap(tmp_path: Path, output_dir: Path):
    """compute_base_level_metrics should treat touching half-open intervals as non-overlapping."""
    predicted_geese = tmp_path / "predicted_touching.geese"
    predicted_geese.write_text(
        "#name\tatom_nr\tclass\tstrand\tstart\tend\n"
        "sequence1\t1\t1\t+\t0\t10\n"
    )

    true_geese = tmp_path / "true_touching.geese"
    true_geese.write_text(
        "#name\tatom_nr\tclass\tstrand\tstart\tend\n"
        "sequence1\t1\t1\t+\t10\t20\n"
    )

    metrics = compute_base_level_metrics(
        predicted_geese=predicted_geese,
        true_geese=true_geese,
        output_directory=output_dir,
        per_class=False
    )

    predicted_status = pd.read_csv(output_dir / "base_predicted_statuses.tsv", sep="\t")
    true_status = pd.read_csv(output_dir / "base_true_statuses.tsv", sep="\t")

    assert metrics == 0.0
    assert predicted_status.to_dict("records") == [
        {
            "name": "sequence1",
            "atom_nr": 1,
            "class": 1,
            "strand": "+",
            "start": 0,
            "end": 10,
            "status": "FP",
        }
    ]
    assert true_status.to_dict("records") == [
        {
            "name": "sequence1",
            "atom_nr": 1,
            "class": 1,
            "strand": "+",
            "start": 10,
            "end": 20,
            "status": "FN",
        }
    ]


# --------------------------------------------------------------------------------------
# Test: basic base-level scoring (overall) with TSV content
# --------------------------------------------------------------------------------------
def test_compute_base_level_metrics(tmp_path: Path, output_dir: Path):
    """compute_base_level_metrics should compute overall base-level metrics and write expected TSV rows."""
    predicted_geese, true_geese = create_partial_overlap_geese(tmp_path)

    metrics = compute_base_level_metrics(
        predicted_geese=predicted_geese,
        true_geese=true_geese,
        output_directory=output_dir,
        per_class=False
    )

    predicted_status = pd.read_csv(output_dir / "base_predicted_statuses.tsv", sep="\t")
    true_status = pd.read_csv(output_dir / "base_true_statuses.tsv", sep="\t")
    overall_metrics = pd.read_csv(output_dir / "base_metrics_overall.tsv", sep="\t")

    assert metrics == 0.5
    assert predicted_status.to_dict("records") == [
        {
            "name": "sequence1",
            "atom_nr": 1,
            "class": 1,
            "strand": "+",
            "start": 0,
            "end": 10,
            "status": "FP",
        },
        {
            "name": "sequence1",
            "atom_nr": 1,
            "class": 1,
            "strand": "+",
            "start": 10,
            "end": 20,
            "status": "TP",
        },
    ]
    assert true_status.to_dict("records") == [
        {
            "name": "sequence1",
            "atom_nr": 1,
            "class": 1,
            "strand": "+",
            "start": 10,
            "end": 20,
            "status": "TP",
        },
        {
            "name": "sequence1",
            "atom_nr": 1,
            "class": 1,
            "strand": "+",
            "start": 20,
            "end": 30,
            "status": "FN",
        },
    ]
    assert overall_metrics.to_dict("records") == [
        {
            "TP": 10,
            "FP": 10,
            "FN": 10,
            "Precision": 0.5,
            "Recall": 0.5,
            "F1-score": 0.5,
        }
    ]


# --------------------------------------------------------------------------------------
# Test: basic base-level scoring (per-class)
# --------------------------------------------------------------------------------------
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
    assert metrics[0]["F1-score"] == 1 / 3
    assert metrics[1]["Class"] == 2
    assert metrics[1]["F1-score"] == 0.0

    predicted_file = output_dir / "base_predicted_statuses.tsv"
    true_file = output_dir / "base_true_statuses.tsv"
    out_file = output_dir / "base_metrics_per_class.tsv"

    assert predicted_file.is_file()
    assert true_file.is_file()
    assert out_file.is_file()


# --------------------------------------------------------------------------------------
# Test: partial class overlap across different classes
# --------------------------------------------------------------------------------------
def test_partial_class_overlap_different_classes(tmp_path: Path, output_dir: Path):
    """
    compute_base_level_metrics should compute metrics when prediction and true values overlap across different classes.
    """
    predicted_geese = tmp_path / "predicted_partial_overlap.geese"
    predicted_geese.write_text(
        "#name\tatom_nr\tclass\tstrand\tstart\tend\n"
        "sequence1\t1\t1\t+\t0\t20\n"
    )

    true_geese = tmp_path / "true_partial_overlap.geese"
    true_geese.write_text(
        "#name\tatom_nr\tclass\tstrand\tstart\tend\n"
        "sequence1\t1\t2\t+\t10\t30\n"
    )

    metrics = compute_base_level_metrics(
        predicted_geese=predicted_geese,
        true_geese=true_geese,
        output_directory=output_dir,
        per_class=False
    )

    assert metrics == 0.0


# --------------------------------------------------------------------------------------
# Test: full class overlap
# --------------------------------------------------------------------------------------
def test_full_class_overlap(tmp_path: Path, output_dir: Path):
    """compute_base_level_metrics should return 1.0 when prediction and true values are identical."""
    predicted_geese = tmp_path / "predicted_full_overlap.geese"
    predicted_geese.write_text(
        "#name\tatom_nr\tclass\tstrand\tstart\tend\n"
        "sequence1\t1\t1\t+\t0\t10\n"
        "sequence1\t2\t1\t+\t10\t20\n"
    )

    true_geese = tmp_path / "true_full_overlap.geese"
    true_geese.write_text(
        "#name\tatom_nr\tclass\tstrand\tstart\tend\n"
        "sequence1\t1\t1\t+\t0\t10\n"
        "sequence1\t2\t1\t+\t10\t20\n"
    )

    metrics = compute_base_level_metrics(
        predicted_geese=predicted_geese,
        true_geese=true_geese,
        output_directory=output_dir,
        per_class=False
    )

    assert metrics == 1.0


# --------------------------------------------------------------------------------------
# Test: empty GEESE files
# --------------------------------------------------------------------------------------
def test_empty_geese(tmp_path: Path, output_dir: Path):
    """compute_base_level_metrics should correctly calculate 0.0 TP, FP, FN."""
    empty_df = pd.DataFrame(columns=pd.Index(["name", "atom_nr", "class", "strand", "start", "end"]))

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
