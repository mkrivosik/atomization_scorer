"""
Tests for the plotting_utils() function.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import pytest

from atomization_scorer.visualization import plotting_utils as pu


# --------------------------------------------------------------------------------------
# Test: normalize_output_format accepts supported formats and normalizes case
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("output_format", "expected"),
    [
        ("png", "png"),
        ("PNG", "png"),
        ("Svg", "svg"),
        ("PDF", "pdf"),
    ],
)
def test_normalize_output_format_supported_values(output_format: str, expected: str):
    """normalize_output_format should lowercase supported figure formats."""
    assert pu.normalize_output_format(output_format) == expected


# --------------------------------------------------------------------------------------
# Test: normalize_output_format rejects unsupported formats
# --------------------------------------------------------------------------------------
def test_normalize_output_format_unsupported_value():
    """normalize_output_format should raise for unsupported figure formats."""
    with pytest.raises(ValueError, match="Unsupported output format 'jpg'"):
        pu.normalize_output_format("jpg")


# --------------------------------------------------------------------------------------
# Test: save_figure saves the figure using the normalized suffix
# --------------------------------------------------------------------------------------
def test_save_figure_writes_normalized_output_file(tmp_path: Path):
    """save_figure should write the figure with a normalized output suffix."""
    fig, _ = plt.subplots()
    output_path = tmp_path / "atomization_plot"

    pu.save_figure(fig=fig, output_path=output_path, output_format="PNG", dpi=150)

    assert output_path.with_suffix(".png").is_file()
    plt.close(fig)


# --------------------------------------------------------------------------------------
# Test: save_figure rejects unsupported formats before saving
# --------------------------------------------------------------------------------------
def test_save_figure_rejects_unsupported_output_format(tmp_path: Path):
    """save_figure should raise ValueError for unsupported output formats."""
    fig, _ = plt.subplots()
    output_path = tmp_path / "atomization_plot"

    with pytest.raises(ValueError, match="Unsupported output format 'jpg'"):
        pu.save_figure(fig=fig, output_path=output_path, output_format="jpg", dpi=150)

    assert not output_path.with_suffix(".jpg").exists()
    plt.close(fig)


# --------------------------------------------------------------------------------------
# Test: get_sorted_intervals filters the selected genome and sorts its intervals
# --------------------------------------------------------------------------------------
def test_get_sorted_intervals_filters_and_sorts_selected_genome():
    """get_sorted_intervals should return sorted intervals only for the requested genome."""
    df = pd.DataFrame(
        {
            "name": ["genome_b", "genome_a", "genome_a", "genome_b"],
            "start": [30, 40, 10, 5],
            "end": [50, 60, 20, 9],
        }
    )

    assert pu.get_sorted_intervals(df, genome_name="genome_a", genome_length=100, label="True") == [
        (10, 20),
        (40, 60),
    ]


# --------------------------------------------------------------------------------------
# Test: get_sorted_intervals returns an empty list when the genome is absent
# --------------------------------------------------------------------------------------
def test_get_sorted_intervals_returns_empty_list_for_missing_genome():
    """get_sorted_intervals should return an empty list when the genome has no intervals."""
    df = pd.DataFrame(
        {
            "name": ["genome_a"],
            "start": [10],
            "end": [20],
        }
    )

    assert pu.get_sorted_intervals(df, genome_name="missing", genome_length=100, label="True") == []


# --------------------------------------------------------------------------------------
# Test: get_sorted_intervals rejects invalid interval coordinates
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("start", "end", "genome_length", "message"),
    [
        (-1, 10, 100, "contains a negative coordinate"),
        (5, -1, 100, "contains a negative coordinate"),
        (10, 10, 100, "must satisfy start < end"),
        (20, 10, 100, "must satisfy start < end"),
        (10, 101, 100, "ends outside genome length 100"),
    ],
)
def test_get_sorted_intervals_rejects_invalid_coordinates(
    start: int,
    end: int,
    genome_length: int,
    message: str,
):
    """get_sorted_intervals should raise ValueError for invalid interval coordinates."""
    df = pd.DataFrame(
        {
            "name": ["genome_a"],
            "start": [start],
            "end": [end],
        }
    )

    with pytest.raises(ValueError, match=message):
        pu.get_sorted_intervals(df, genome_name="genome_a", genome_length=genome_length, label="True")


# --------------------------------------------------------------------------------------
# Test: get_sorted_intervals warns on overlaps but still returns sorted intervals
# --------------------------------------------------------------------------------------
def test_get_sorted_intervals_warns_on_overlap(caplog: pytest.LogCaptureFixture):
    """get_sorted_intervals should log a warning when sorted intervals overlap."""
    df = pd.DataFrame(
        {
            "name": ["genome_a", "genome_a"],
            "start": [20, 10],
            "end": [30, 25],
        }
    )

    with caplog.at_level("WARNING"):
        intervals = pu.get_sorted_intervals(
            df,
            genome_name="genome_a",
            genome_length=100,
            label="Predicted",
        )

    assert intervals == [(10, 25), (20, 30)]
    assert "Predicted intervals for genome 'genome_a' overlap near (20, 30)" in caplog.text


# --------------------------------------------------------------------------------------
# Test: split_interval_for_rows preserves wrapped fragments and true boundary flags
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("start", "end", "line_length", "expected_fragments"),
    [
        (
            1000,
            2000,
            5000,
            [(0, 1000, 2000, True, True)],
        ),
        (
            3000,
            13000,
            5000,
            [
                (0, 3000, 5000, True, False),
                (1, 0, 5000, False, False),
                (2, 0, 3000, False, True),
            ],
        ),
        (
            5000,
            7000,
            5000,
            [(1, 0, 2000, True, True)],
        ),
        (
            2000,
            5000,
            5000,
            [(0, 2000, 5000, True, True)],
        ),
    ],
)
def test_split_interval_for_rows_expected_fragments(
    start: int,
    end: int,
    line_length: int,
    expected_fragments: list[tuple[int, int, int, bool, bool]],
):
    """split_interval_for_rows should yield row-local fragments with correct boundary flags."""
    # noinspection PyTypeChecker
    assert list(pu.split_interval_for_rows(start, end, line_length)) == expected_fragments


# --------------------------------------------------------------------------------------
# Test: compute_gap_segments returns expected row-local gaps
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("intervals", "row", "line_length", "genome_length", "expected_gaps"),
    [
        ([], 0, 5000, 5000, [(0, 5000)]),
        ([(0, 5000)], 0, 5000, 5000, []),
        ([(1000, 3000)], 0, 5000, 5000, [(0, 1000), (3000, 5000)]),
        ([(4500, 6500)], 1, 5000, 10000, [(1500, 5000)]),
        ([(5200, 7000), (7600, 8300)], 1, 5000, 10000, [(0, 200), (2000, 2600), (3300, 5000)]),
        ([(0, 4900), (5400, 9999), (10000, 18000)], 1, 5000, 20000, [(0, 400), (4999, 5000)]),
        ([(0, 4900), (5400, 9999), (10000, 18000)], 3, 5000, 20000, [(3000, 5000)]),
    ],
)
def test_compute_gap_segments_expected_gaps(
    intervals: list[tuple[int, int]],
    row: int,
    line_length: int,
    genome_length: int,
    expected_gaps: list[tuple[int, int]],
):
    """compute_gap_segments should return uncovered row-local segments for wrapped rows."""
    assert pu.compute_gap_segments(intervals, row, line_length, genome_length) == expected_gaps
