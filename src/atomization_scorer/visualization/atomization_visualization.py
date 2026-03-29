"""
atomization_visualization.py

Wrapped atomization visualization of predicted versus true atomization.

Functions
---------
_draw_baseline      : Draw uncovered genome baseline segments for a wrapped track.
_draw_atom_track    : Draw colored atom segments plus true start/end boundary markers.
plot_atomization    : Generate one wrapped atomization figure per genome.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
from __future__ import annotations
import logging
import math
from pathlib import Path

import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.collections import LineCollection, PatchCollection

from atomization_scorer import read_fasta, read_geese

from .plotting_utils import (
    compute_gap_segments,
    get_sorted_intervals,
    normalize_output_format,
    save_figure,
    split_interval_for_rows,
)

# --------------------------------------------------------------------------------------
# Logging
# --------------------------------------------------------------------------------------
logger = logging.getLogger(__name__)


# --------------------------------------------------------------------------------------
# Helper: draw uncovered genome baseline for one wrapped track
# --------------------------------------------------------------------------------------
def _draw_baseline(
    ax: Axes,
    intervals: list[tuple[int, int]],
    genome_length: int,
    line_length: int,
    track_y: float,
) -> None:
    """
    Draw black baseline segments for genome regions not covered by atom intervals.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis where the baseline segments are drawn.
    intervals : list[tuple[int, int]]
        Atom intervals for one track.
    genome_length : int
        Total length of the genome.
    line_length : int
        Number of genome bases shown in one wrapped row.
    track_y : float
        Base y-position of the track inside each wrapped row.

    Returns
    -------
    None
    """
    lines = []
    n_rows = math.ceil(genome_length / line_length)
    for row in range(n_rows):
        y = track_y + row * 3.0
        for gap_start, gap_end in compute_gap_segments(intervals, row, line_length, genome_length):
            lines.append([[gap_start, y], [gap_end, y]])

    if lines:
        ax.add_collection(
            LineCollection(
                lines,
                colors="black",
                linewidths=2.0,
                zorder=1,
            )
        )


# --------------------------------------------------------------------------------------
# Helper: draw atom segments and boundary markers for one wrapped track
# --------------------------------------------------------------------------------------
def _draw_atom_track(
    ax: Axes,
    intervals: list[tuple[int, int]],
    line_length: int,
    track_y: float,
    color: str,
) -> None:
    """
    Draw colored atom segments and their true start and end boundary markers.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis where the atom segments and boundary markers are drawn.
    intervals : list[tuple[int, int]]
        Atom intervals for one track.
    line_length : int
        Number of genome bases shown in one wrapped row.
    track_y : float
        Base y-position of the track inside each wrapped row.
    color : str
        Color used for the atom segments and boundary markers.

    Returns
    -------
    None
    """
    rectangles = []
    boundaries = []

    for start, end in intervals:
        # noinspection PyTypeChecker
        for fragment in split_interval_for_rows(start, end, line_length):
            row, x_start, x_end, is_true_start, is_true_end = fragment
            y = track_y + row * 3.0
            rectangles.append(
                patches.Rectangle(
                    (x_start, y - 0.075),
                    x_end - x_start,
                    0.15,
                )
            )
            if is_true_start:
                boundaries.append([[x_start, y - 0.175], [x_start, y + 0.175]])
            if is_true_end:
                boundaries.append([[x_end, y - 0.175], [x_end, y + 0.175]])

    if rectangles:
        ax.add_collection(
            PatchCollection(
                rectangles,
                facecolor=color,
                edgecolor=color,
                alpha=1.0,
                zorder=3
            )
        )
    if boundaries:
        ax.add_collection(
            LineCollection(
                boundaries,
                colors=color,
                linewidths=2.0,
                zorder=4,
                clip_on=False,
                capstyle="projecting",
            )
        )


# --------------------------------------------------------------------------------------
# Atomization Visualization
# --------------------------------------------------------------------------------------
def plot_atomization(
    genomes_file: Path,
    true_atoms_file: Path,
    predicted_atoms_file: Path,
    output_directory: Path,
    figure_width: float = 12.0,
    dpi: int = 150,
    target_rows: int = 20,
    min_bases_per_row: int = 10_000,
    max_bases_per_row: int = 250_000,
    true_color: str = "#2C7FB8",
    predicted_color: str = "#F28E2B",
    output_format: str = "png",
) -> None:
    """
    Generate one wrapped atomization figure per genome comparing predicted and true atomization.

    Parameters
    ----------
    genomes_file : Path
        Input genomes FASTA file containing the genome sequences.
    true_atoms_file : Path
        Input GEESE file containing the true atomization.
    predicted_atoms_file : Path
        Input GEESE file containing the predicted atomization.
    output_directory : Path
        Path to the output directory where figures are stored.
    figure_width : float, optional, default=12.0
        Width of the output figure in inches.
    dpi : int, optional, default=150
        Resolution of the saved figure.
    target_rows : int, optional, default=20
        Target number of wrapped rows per genome.
    min_bases_per_row : int, optional, default=10_000
        Minimum number of genome bases shown in one wrapped row.
    max_bases_per_row : int, optional, default=250_000
        Maximum number of genome bases shown in one wrapped row.
    true_color : str, optional, default="#2C7FB8"
        Color used for the true atomization track.
    predicted_color : str, optional, default="#F28E2B"
        Color used for the predicted atomization track.
    output_format : str, optional, default="png"
        Output figure format.

    Raises
    ------
    FileNotFoundError
        Raised if genomes_file, true_atoms_file, or predicted_atoms_file do not exist.
    ValueError
        Raised if target_rows, min_bases_per_row, or max_bases_per_row are invalid.

    Returns
    -------
    None
    """
    if not genomes_file.is_file():
        raise FileNotFoundError(f"Genomes FASTA file not found: {genomes_file}")
    if not true_atoms_file.is_file():
        raise FileNotFoundError(f"True atoms file not found: {true_atoms_file}")
    if not predicted_atoms_file.is_file():
        raise FileNotFoundError(f"Predicted atoms file not found: {predicted_atoms_file}")
    if target_rows <= 0:
        raise ValueError("target_rows must be a positive integer.")
    if min_bases_per_row <= 0:
        raise ValueError("min_bases_per_row must be a positive integer.")
    if max_bases_per_row <= 0:
        raise ValueError("max_bases_per_row must be a positive integer.")
    if min_bases_per_row > max_bases_per_row:
        raise ValueError("min_bases_per_row must be less than or equal to max_bases_per_row.")

    normalized_format = normalize_output_format(output_format)
    output_directory.mkdir(parents=True, exist_ok=True)

    logger.info(
        (
            "Generating wrapped atomization visualizations from "
            "genomes=%s true_atoms=%s predicted_atoms=%s into %s as %s"
        ),
        genomes_file,
        true_atoms_file,
        predicted_atoms_file,
        output_directory,
        normalized_format,
    )

    genome_dictionary = read_fasta(genomes_file)
    df_true = read_geese(true_atoms_file)
    df_predicted = read_geese(predicted_atoms_file)

    for genome_name, sequence in genome_dictionary.items():
        genome_length = len(sequence)
        line_length = math.ceil(genome_length / target_rows)
        line_length = max(min_bases_per_row, line_length)
        line_length = min(max_bases_per_row, line_length, genome_length)
        n_rows = max(1, math.ceil(genome_length / line_length))

        true_intervals = get_sorted_intervals(
            df_true,
            genome_name,
            genome_length,
            "True",
        )
        predicted_intervals = get_sorted_intervals(
            df_predicted,
            genome_name,
            genome_length,
            "Predicted",
        )

        logger.info(
            (
                "Rendering wrapped atomization visualization for genome=%s "
                "length=%s true_intervals=%s predicted_intervals=%s rows=%s "
                "line_length=%s target_rows=%s min_bases_per_row=%s "
                "max_bases_per_row=%s"
            ),
            genome_name,
            genome_length,
            len(true_intervals),
            len(predicted_intervals),
            n_rows,
            line_length,
            target_rows,
            min_bases_per_row,
            max_bases_per_row,
        )

        figure_height = max(3.0, 1.5 + n_rows * 1.8)
        fig, ax = plt.subplots(figsize=(figure_width, figure_height))

        _draw_baseline(ax, predicted_intervals, genome_length, line_length, 1.0)
        _draw_baseline(ax, true_intervals, genome_length, line_length, 2.0)
        _draw_atom_track(ax, predicted_intervals, line_length, 1.0, predicted_color)
        _draw_atom_track(ax, true_intervals, line_length, 2.0, true_color)

        ax.set_xlim(-0.5, line_length + 0.5)
        ax.set_ylim(2.50 + (n_rows - 1) * 3.0, 0.5)

        ticks = []
        labels = []
        for row in range(n_rows):
            row_start = row * line_length
            row_end = min(row_start + line_length, genome_length) - 1
            ticks.append(1.0 + row * 3.0)
            labels.append(f"Predicted {row_start}-{row_end}")
            ticks.append(2.0 + row * 3.0)
            labels.append(f"True {row_start}-{row_end}")

        ax.set_yticks(ticks)
        ax.set_yticklabels(labels)
        ax.set_xlabel("Genome position within wrapped row")
        ax.set_title(f"Atomization: {genome_name}")

        fig.tight_layout()
        output_path = output_directory / genome_name
        save_figure(fig, output_path, normalized_format, dpi)
        plt.close(fig)

        logger.info(
            "Saved wrapped atomization visualization for genome=%s to %s.%s",
            genome_name,
            output_path,
            normalized_format,
        )
