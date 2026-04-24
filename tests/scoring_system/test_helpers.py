"""
Tests for the helper.py functions.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
from pathlib import Path

import pandas as pd
import pytest

from atomization_scorer.scoring_system.helpers import (
    _compute_and_write_metrics,  # noqa
    _compute_metrics,  # noqa
    _create_new_row,  # noqa
    _interval_overlap,  # noqa
    _write_metrics_tsv,  # noqa
)


# --------------------------------------------------------------------------------------
# Tests for _compute_and_write_metrics
# --------------------------------------------------------------------------------------
def test_compute_and_write_metrics_overall_basic(output_dir: Path):
    """_compute_and_write_metrics should be computed and written correctly."""
    tp = {1: 10, 2: 5}
    fp = {1: 2}
    fn = {2: 3}

    output_file = output_dir / 'metrics.tsv'

    f1 = _compute_and_write_metrics(
        tp=tp, fp=fp, fn=fn,
        output_file=output_file,
        per_class=False
    )

    assert isinstance(f1, float)
    assert output_file.is_file()

    df = pd.read_csv(output_file, sep="\t")
    assert list(df.columns) == ["TP", "FP", "FN", "Precision", "Recall", "F1-score"]
    assert df.iloc[0]["TP"] == 15
    assert df.iloc[0]["FP"] == 2
    assert df.iloc[0]["FN"] == 3
    assert df.iloc[0]["Precision"] == 15 / 17
    assert df.iloc[0]["Recall"] == 15 / 18
    assert df.iloc[0]["F1-score"] == 6 / 7


def test_compute_and_write_metrics_overall_all_zero(output_dir: Path):
    """_compute_and_write_metrics should return zero scores if all counts are zero."""
    tp = {}
    fp = {}
    fn = {}

    output_file = output_dir / "overall_zero.tsv"

    f1 = _compute_and_write_metrics(
        tp=tp, fp=fp, fn=fn,
        output_file=output_file,
        per_class=False
    )

    assert f1 == 0.0

    df = pd.read_csv(output_file, sep="\t")
    assert df.iloc[0]["Precision"] == 0.0
    assert df.iloc[0]["Recall"] == 0.0
    assert df.iloc[0]["F1-score"] == 0.0


def test_compute_and_write_metrics_per_class_basic(output_dir: Path):
    """_compute_and_write_metrics should be computed for all present classes."""
    tp = {1: 10}
    fp = {1: 2, 2: 3}
    fn = {2: 5}

    output_file = output_dir / "per_class.tsv"

    result = _compute_and_write_metrics(
        tp=tp, fp=fp, fn=fn,
        output_file=output_file,
        per_class=True
    )

    assert isinstance(result, list)
    assert output_file.is_file()
    assert len(result) == 2
    assert result == [
        {"Class": 1, "F1-score": 10 / 11},
        {"Class": 2, "F1-score": 0.0},
    ]

    df = pd.read_csv(output_file, sep="\t")
    assert df.to_dict("records") == [
        {
            "Class": 1,
            "TP": 10,
            "FP": 2,
            "FN": 0,
            "Precision": pytest.approx(10 / 12),
            "Recall": pytest.approx(1.0),
            "F1-score": pytest.approx(10 / 11),
        },
        {
            "Class": 2,
            "TP": 0,
            "FP": 3,
            "FN": 5,
            "Precision": pytest.approx(0.0),
            "Recall": pytest.approx(0.0),
            "F1-score": pytest.approx(0.0),
        },
    ]


def test_compute_and_write_metrics_classes_sorted(output_dir: Path):
    """_compute_and_write_metrics: per-class metrics should be sorted by atomization class."""
    tp = {10: 1, 2: 1}
    fp = {}
    fn = {}

    output_file = output_dir / "sorted.tsv"

    _compute_and_write_metrics(
        tp=tp, fp=fp, fn=fn,
        output_file=output_file,
        per_class=True
    )

    df = pd.read_csv(output_file, sep="\t")
    assert list(df["Class"]) == [2, 10]


def test_compute_and_write_metrics_return_format(output_dir: Path):
    """_compute_and_write_metrics: return value should only contain Class and F1-score."""
    tp = {1: 3}
    fp = {}
    fn = {}

    output_file = output_dir / "return_format.tsv"

    result = _compute_and_write_metrics(
        tp=tp, fp=fp, fn=fn,
        output_file=output_file,
        per_class=True
    )

    assert result == [{"Class": 1, "F1-score": 1.0}]


def test_compute_and_write_metrics_per_class_fn_only_class(output_dir: Path):
    """_compute_and_write_metrics should include classes that appear only in FN."""
    tp = {}
    fp = {}
    fn = {7: 4}

    output_file = output_dir / "fn_only.tsv"

    result = _compute_and_write_metrics(
        tp=tp, fp=fp, fn=fn,
        output_file=output_file,
        per_class=True
    )

    df = pd.read_csv(output_file, sep="\t")

    assert result == [{"Class": 7, "F1-score": 0.0}]
    assert df.to_dict("records") == [
        {
            "Class": 7,
            "TP": 0,
            "FP": 0,
            "FN": 4,
            "Precision": 0.0,
            "Recall": 0.0,
            "F1-score": 0.0,
        }
    ]


# --------------------------------------------------------------------------------------
# Tests for _compute_metrics
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("tp", "fp", "fn", "expected_precision", "expected_recall", "expected_f1"),
    [
        (5, 5, 5, 0.5, 0.5, 0.5),
        (10, 0, 0, 1.0, 1.0, 1.0),
        (0, 0, 0, 0.0, 0.0, 0.0),
        (0, 5, 0, 0.0, 0.0, 0.0),
        (0, 0, 5, 0.0, 0.0, 0.0),
    ],
    ids=[
        "mixed-values",
        "perfect-prediction",
        "all-zero",
        "false-positives-only",
        "false-negatives-only",
    ],
)
def test_compute_metrics_cases(
    tp: int,
    fp: int,
    fn: int,
    expected_precision: float,
    expected_recall: float,
    expected_f1: float,
):
    """_compute_metrics should compute the expected values for simple scoring cases."""
    precision, recall, f1 = _compute_metrics(tp=tp, fp=fp, fn=fn)

    assert precision == expected_precision
    assert recall == expected_recall
    assert f1 == expected_f1


@pytest.mark.parametrize(
    ("tp", "fp", "fn", "expected_precision", "expected_recall", "expected_f1"),
    [
        (4, 0, 2, 1.0, 4 / 6, 4 / 5),
        (4, 2, 0, 4 / 6, 1.0, 4 / 5),
    ],
    ids=[
        "perfect-precision-partial-recall",
        "partial-precision-perfect-recall",
    ],
)
def test_compute_metrics_asymmetric_cases(
    tp: int,
    fp: int,
    fn: int,
    expected_precision: float,
    expected_recall: float,
    expected_f1: float,
):
    """_compute_metrics should handle one-sided partial precision and recall correctly."""
    precision, recall, f1 = _compute_metrics(tp=tp, fp=fp, fn=fn)

    assert precision == expected_precision
    assert recall == expected_recall
    assert f1 == expected_f1


# --------------------------------------------------------------------------------------
# Tests for _write_metrics_tsv
# --------------------------------------------------------------------------------------
def test_write_metrics_tsv_creates_file(output_dir: Path):
    """_write_metrics_tsv should create a TSV file."""
    df = pd.DataFrame(
        [[1, 2, 3]],
        columns=pd.Index(["TP", "FP", "FN"])
    )

    output_file = output_dir / "metrics.tsv"

    _write_metrics_tsv(df=df, output_path=output_file)

    assert output_file.is_file()


def test_write_metrics_tsv_creates_parent_directory(tmp_path: Path):
    """_write_metrics_tsv should create parent directories if they do not exist."""
    df = pd.DataFrame(
        [[1, 2]],
        columns=pd.Index(["Precision", "Recall"])
    )

    output_file = tmp_path / "nested" / "directory" / "metrics.tsv"

    _write_metrics_tsv(df=df, output_path=output_file)

    assert output_file.is_file()
    assert output_file.parent.exists()


def test_write_metrics_tsv_empty_dataframe(output_dir: Path):
    """_write_metrics_tsv should correctly write an empty DataFrame."""
    df = pd.DataFrame(columns=pd.Index(["A", "B", "C"]))

    output_file = output_dir / "empty.tsv"
    _write_metrics_tsv(df=df, output_path=output_file)

    assert output_file.is_file()

    read_df = pd.read_csv(output_file, sep="\t")
    assert read_df.empty
    assert list(read_df.columns) == ["A", "B", "C"]


def test_write_metrics_tsv_overwrites_existing_file(output_dir: Path):
    """_write_metrics_tsv should overwrite an existing file."""
    output_file = output_dir / "overwrite.tsv"

    df1 = pd.DataFrame([[1]], columns=pd.Index(["A"]))
    df2 = pd.DataFrame([[2]], columns=pd.Index(["A"]))

    _write_metrics_tsv(df=df1, output_path=output_file)
    _write_metrics_tsv(df=df2, output_path=output_file)

    read_df = pd.read_csv(output_file, sep="\t")
    assert read_df.iloc[0]["A"] == 2


# --------------------------------------------------------------------------------------
# Tests for _create_new_row
# --------------------------------------------------------------------------------------
def test_create_new_row_preserves_metadata_and_replaces_coordinates():
    """_create_new_row should keep metadata while replacing the original coordinates."""
    row = pd.Series({
        "name": "sequence_2",
        "atom_nr": 9,
        "class": 4,
        "strand": "-",
        "start": 100,
        "end": 250
    })

    new_row = _create_new_row(
        row=row,
        start=120,
        end=140,
        status="FN"
    )

    assert new_row == {
        "name": "sequence_2",
        "atom_nr": 9,
        "class": 4,
        "strand": "-",
        "start": 120,
        "end": 140,
        "status": "FN",
    }


# --------------------------------------------------------------------------------------
# Tests for _interval_overlap
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("start1", "end1", "start2", "end2", "expected_overlap"),
    [
        (0, 10, 20, 30, 0.0),
        (5, 15, 5, 15, 1.0),
        (0, 10, 6, 15, 4 / 15),
        (0, 20, 5, 10, 5 / 20),
        (0, 10, 10, 20, 0.0),
    ],
    ids=[
        "no-overlap",
        "identical-intervals",
        "partial-overlap",
        "interval-inside-another",
        "touching-boundary",
    ],
)
def test_interval_overlap_cases(
    start1: int,
    end1: int,
    start2: int,
    end2: int,
    expected_overlap: float,
):
    """_interval_overlap should compute the expected overlap ratio for common interval cases."""
    overlap = _interval_overlap(start1=start1, end1=end1, start2=start2, end2=end2)
    assert overlap == expected_overlap
