"""
plotting_utils.py

Shared helper functions for atomization visualization.

Functions
---------
normalize_output_format : Validate and normalize a requested figure output format.
save_figure             : Save a matplotlib figure using the requested format.
get_sorted_intervals    : Extract, sort, and validate intervals for one genome.
split_interval_for_rows : Split an interval across wrapped rows while preserving boundaries.
compute_gap_segments    : Compute uncovered baseline segments for one wrapped row.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
from __future__ import annotations
from pathlib import Path
from typing import Iterable

import pandas as pd
import logging
from matplotlib.figure import Figure

# --------------------------------------------------------------------------------------
# Logging
# --------------------------------------------------------------------------------------
logger = logging.getLogger(__name__)

# --------------------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------------------
SUPPORTED_OUTPUT_FORMATS = {"png", "svg", "pdf"}


# --------------------------------------------------------------------------------------
# Plotting Utils Functions
# --------------------------------------------------------------------------------------
def normalize_output_format(output_format: str) -> str:
    """
    Validate and normalize a requested figure output format.

    Parameters
    ----------
    output_format : str
        Requested figure output format.

    Raises
    ------
    ValueError
        Raised if output_format is not one of the supported figure formats.

    Returns
    -------
    str
        Normalized lowercase output format.
    """
    normalized = output_format.lower()
    if normalized not in SUPPORTED_OUTPUT_FORMATS:
        supported = ", ".join(sorted(SUPPORTED_OUTPUT_FORMATS))
        raise ValueError(
            f"Unsupported output format '{output_format}'. Supported formats: {supported}."
        )
    return normalized


def save_figure(fig: Figure, output_path: Path, output_format: str, dpi: int) -> None:
    """
    Save a matplotlib figure using the requested format.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        A matplotlib figure to save.
    output_path : Path
        Output file path without the figure format suffix.
    output_format : str
        Requested figure output format.
    dpi : int
        Resolution used when saving the figure.

    Raises
    ------
    ValueError
        Raised if output_format is not one of the supported figure formats.

    Returns
    -------
    None
    """
    normalized_format = normalize_output_format(output_format)
    fig.savefig(
        output_path.with_suffix(f".{normalized_format}"),
        format=normalized_format,
        dpi=dpi,
        bbox_inches="tight",
    )


def get_sorted_intervals(
    df: pd.DataFrame,
    genome_name: str,
    genome_length: int,
    label: str,
) -> list[tuple[int, int]]:
    """
    Extract, sort, and validate intervals for a single genome.

    Parameters
    ----------
    df : pandas.DataFrame
        Table containing atom intervals for one or more genomes.
    genome_name : str
        Name of the genome whose intervals are extracted.
    genome_length : int
        Total length of the genome.
    label : str
        Label used in warning and error messages to identify the interval source.

    Raises
    ------
    ValueError
        Raised if an interval contains a negative coordinate, does not satisfy
        start < end, or ends outside the genome length.

    Returns
    -------
    list[tuple[int, int]]
        Sorted list of validated (start, end) atom intervals for the genome.
    """
    intervals = [
        (int(start), int(end))
        for start, end in zip(
            df.loc[df["name"] == genome_name, "start"],
            df.loc[df["name"] == genome_name, "end"],
        )
    ]
    intervals.sort()

    previous_end = None
    for start, end in intervals:
        if start < 0 or end < 0:
            raise ValueError(f"{label} interval for genome '{genome_name}' contains a negative coordinate.")
        if start >= end:
            raise ValueError(
                f"{label} interval for genome '{genome_name}' must satisfy start < end, got ({start}, {end})."
            )
        if end > genome_length:
            raise ValueError(
                f"{label} interval for genome '{genome_name}' ends outside genome length {genome_length}: ({start}, {end})."
            )
        if previous_end is not None and start < previous_end:
            logger.warning(
                "%s intervals for genome '%s' overlap near (%s, %s); plotting will continue.",
                label,
                genome_name,
                start,
                end,
            )
        previous_end = end

    return intervals


def split_interval_for_rows(
    start: int,
    end: int,
    line_length: int,
) -> Iterable[tuple[int, int, int, bool, bool]]:
    """
    Split a single interval across wrapped rows while preserving true boundaries.

    Parameters
    ----------
    start : int
        Start coordinate of the interval.
    end : int
        End coordinate of the interval.
    line_length : int
        Number of genome bases shown in one wrapped row.

    Yields
    ------
    tuple[int, int, int, bool, bool]
        Tuples of the form (row_index, x_start, x_end, is_true_start, is_true_end),
        where x_start and x_end are row-local coordinates and the boolean values
        indicate whether the fragment contains the true interval start or end.

    Returns
    -------
    Iterable[tuple[int, int, int, bool, bool]]
        Iterable of wrapped interval fragments for plotting.
    """
    current = start
    while current < end:
        row_index = current // line_length
        row_start = row_index * line_length
        row_end = row_start + line_length
        fragment_end = min(end, row_end)
        yield (
            row_index,
            current - row_start,
            fragment_end - row_start,
            current == start,
            fragment_end == end,
        )
        current = fragment_end


def compute_gap_segments(
    intervals: list[tuple[int, int]],
    row: int,
    line_length: int,
    genome_length: int,
) -> list[tuple[int, int]]:
    """
    Compute uncovered baseline segments for one wrapped row.

    Parameters
    ----------
    intervals : list[tuple[int, int]]
        Sorted atom intervals in genome coordinates.
    row : int
        Zero-based wrapped row index.
    line_length : int
        Number of genome bases shown in one wrapped row.
    genome_length : int
        Total length of the genome.

    Returns
    -------
    list[tuple[int, int]]
        Uncovered segments in row-local coordinates for the selected wrapped row.
    """
    row_start = row * line_length
    row_end = min(row_start + line_length, genome_length)
    intervals_in_row = []

    for start, end in intervals:
        if end <= row_start or start >= row_end:
            continue
        intervals_in_row.append((max(start, row_start), min(end, row_end)))

    current = row_start
    gaps = []
    for start, end in intervals_in_row:
        if start > current:
            gaps.append((current - row_start, start - row_start))
        current = max(current, end)

    if current < row_end:
        gaps.append((current - row_start, row_end - row_start))

    return gaps
