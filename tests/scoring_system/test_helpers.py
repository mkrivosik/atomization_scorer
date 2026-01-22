"""
Tests for the helper.py functions.
"""

from pathlib import Path

import pandas as pd

from atomization_scorer.scoring_system.helpers import (
    _compute_and_write_metrics,  # noqa
    _compute_metrics,  # noqa
    _create_new_row,  # noqa
    _interval_overlap,  # noqa
    _write_metrics_tsv,  # noqa
)


# ---------------------------------------------------------------------
# Tests for _compute_and_write_metrics
# ---------------------------------------------------------------------
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

    df = pd.read_csv(output_file, sep="\t")
    assert set(df["Class"]) == {1, 2}


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


# ---------------------------------------------------------------------
# Tests for _compute_metrics
# ---------------------------------------------------------------------
def test_compute_metrics_mixed_values():
    """_compute_metrics should compute correct values for mixed TP, FP, FN."""
    precision, recall, f1 = _compute_metrics(tp=5, fp=5, fn=5)

    assert precision == 0.5
    assert recall == 0.5
    assert f1 == 0.5


def test_compute_metrics_perfect_prediction():
    """_compute_metrics should return 1.0 for perfect prediction."""
    precision, recall, f1 = _compute_metrics(tp=10, fp=0, fn=0)

    assert precision == 1.0
    assert recall == 1.0
    assert f1 == 1.0


def test_compute_metrics_zero_everything():
    """_compute_metrics should return all zeros when tp, fp, fn are zero."""
    precision, recall, f1 = _compute_metrics(tp=0, fp=0, fn=0)

    assert precision == 0.0
    assert recall == 0.0
    assert f1 == 0.0


def test_compute_metrics_only_false_positives():
    """_compute_metrics should return zero precision, recall, and f1 if only FP exist."""
    precision, recall, f1 = _compute_metrics(tp=0, fp=5, fn=0)

    assert precision == 0.0
    assert recall == 0.0
    assert f1 == 0.0


# ---------------------------------------------------------------------
# Tests for _write_metrics_tsv
# ---------------------------------------------------------------------
def test_write_metrics_tsv_creates_file(output_dir: Path):
    """_write_metrics_tsv should create a TSV file."""
    df = pd.DataFrame(
        [[1, 2, 3]],
        columns=["TP", "FP", "FN"]
    )

    output_file = output_dir / "metrics.tsv"

    _write_metrics_tsv(df=df, output_path=output_file)

    assert output_file.is_file()


def test_write_metrics_tsv_creates_parent_directory(tmp_path: Path):
    """_write_metrics_tsv should create parent directories if they do not exist."""
    df = pd.DataFrame(
        [[1, 2]],
        columns=["Precision", "Recall"]
    )

    output_file = tmp_path / "nested" / "directory" / "metrics.tsv"

    _write_metrics_tsv(df=df, output_path=output_file)

    assert output_file.is_file()
    assert output_file.parent.exists()


def test_write_metrics_tsv_empty_dataframe(output_dir: Path):
    """_write_metrics_tsv should correctly write an empty DataFrame."""
    df = pd.DataFrame(columns=["A", "B", "C"])

    output_file = output_dir / "empty.tsv"
    _write_metrics_tsv(df=df, output_path=output_file)

    assert output_file.is_file()

    read_df = pd.read_csv(output_file, sep="\t")
    assert read_df.empty
    assert list(read_df.columns) == ["A", "B", "C"]


def test_write_metrics_tsv_overwrites_existing_file(output_dir: Path):
    """_write_metrics_tsv should overwrite existing file."""
    output_file = output_dir / "overwrite.tsv"

    df1 = pd.DataFrame([[1]], columns=["A"])
    df2 = pd.DataFrame([[2]], columns=["A"])

    _write_metrics_tsv(df=df1, output_path=output_file)
    _write_metrics_tsv(df=df2, output_path=output_file)

    read_df = pd.read_csv(output_file, sep="\t")
    assert read_df.iloc[0]["A"] == 2


# ---------------------------------------------------------------------
# Tests for _create_new_row
# ---------------------------------------------------------------------
def test_create_new_row_basic():
    """_create_new_row should correctly construct a new interval row."""
    row = pd.Series({
        "name": "sequence",
        "atom_nr": 5,
        "class": 2,
        "strand": "+",
        "start": 0,
        "end": 100
    })

    new_row = _create_new_row(
        row=row,
        start=10,
        end=20,
        status="TP"
    )

    assert isinstance(new_row, dict)
    assert new_row["name"] == "sequence"
    assert new_row["atom_nr"] == 5
    assert new_row["class"] == 2
    assert new_row["strand"] == "+"
    assert new_row["start"] == 10
    assert new_row["end"] == 20
    assert new_row["status"] == "TP"


# ---------------------------------------------------------------------
# Tests for _interval_overlap
# ---------------------------------------------------------------------
def test_interval_overlap_no_overlap():
    """_interval_overlap should return 0.0 if intervals do not overlap."""
    # overlap = 0...10 => 11 bp
    # union = 20...30 => 11 bp
    overlap = _interval_overlap(start1=0, end1=10, start2=20, end2=30)
    assert overlap == 0.0


def test_interval_overlap_identical_intervals():
    """_interval_overlap should return 1.0 for identical intervals."""
    # overlap = 5...15 => 11 bp
    # union = 5...15 => 11 bp
    overlap = _interval_overlap(start1=5, end1=15, start2=5, end2=15)
    assert overlap == 1.0


def test_interval_overlap_partial_overlap():
    """_interval_overlap should correctly compute partial overlap."""
    # overlap = 6...10 => 5 bp
    # union = 0...15 => 16 bp
    overlap = _interval_overlap(start1=0, end1=10, start2=6, end2=15)
    assert overlap == 5 / 16


def test_interval_overlap_interval_inside_another():
    """_interval_overlap should handle one interval fully inside another."""
    # overlap = 5...10 => 6 bp
    # union = 0...20 => 21 bp
    overlap = _interval_overlap(start1=0, end1=20, start2=5, end2=10)
    assert overlap == 6 / 21


def test_interval_overlap_single_base_overlap():
    """_interval_overlap should handle single-base overlap."""
    # overlap at position 10
    # union = 0...20 => 21 bp
    overlap = _interval_overlap(start1=0, end1=10, start2=10, end2=20)
    assert overlap == 1 / 21


def test_interval_overlap_single_base_intervals():
    """_interval_overlap should handle single-base intervals."""
    # overlap at position 5
    # union at position 5
    overlap = _interval_overlap(start1=5, end1=5, start2=5, end2=5)
    assert overlap == 1.0


def test_interval_overlap_adjacent_intervals():
    """_interval_overlap should return 0.0 for adjacent but non-overlapping intervals."""
    overlap = _interval_overlap(start1=0, end1=9, start2=10, end2=20)
    # overlap = 0...9 => 10 bp
    # union = 10...20 => 11 bp
    assert overlap == 0.0
