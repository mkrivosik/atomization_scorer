"""
atomization_visualization.py

Provides functionality to visualize genome atomization for each genome separately.

Modules
-------
plot_genome_atomization : Generates PNG visualizations of true vs predicted atoms per genome.
"""

# ---------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection, LineCollection

from atomization_scorer import read_fasta, read_geese

# ---------------------------------------------------------------------
# Visualization style constants
# ---------------------------------------------------------------------
RECT_HEIGHT = 0.15
VLINE_HEIGHT = 0.35
ROW_SPACING = 3
TARGET_ROWS = 20

# ---------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------

def _wrapped_intervals(
    intervals: list[tuple[int, int]],
    line_length: int,
) -> list[tuple[int, int]]:
    """
    Split genome atom intervals into sub-intervals that fit within a single row.

    This function ensures that any interval spanning multiple rows is divided so that
    each piece fits within `line_length`. This is used to layout atoms across multiple
    rows for visualization.

    Parameters
    ----------
    intervals : list of tuple[int, int]
        List of start-end positions of atoms for the genome.
    line_length : int
        Maximum number of bases per row.

    Returns
    -------
    list of tuple[int, int]
        List of intervals split to fit within row boundaries.
    """
    new_intervals: list[tuple[int, int]] = []
    for start, end in intervals:
        s = start
        while s < end:
            row_start = (s // line_length) * line_length
            row_end = row_start + line_length
            new_end = min(end, row_end)
            new_intervals.append((s, new_end))
            s = new_end
    return new_intervals


def _draw_atoms(
    ax,
    intervals: list[tuple[int, int]],
    line_length: int,
    base_y: float,
    color: str,
) -> None:
    """
    Draw rectangles representing genome atoms along a track, with vertical
    boundary lines marking true atom starts and ends.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The matplotlib axes to draw on.
    intervals : list of tuple[int, int]
        List of start-end positions of atoms for the genome.
    line_length : int
        Number of bases per row (used for wrapping).
    base_y : float
        Y-coordinate of the base track for this set of atoms.
    color : str
        Color of the rectangles and vertical lines.

    Returns
    -------
    None
    """
    rects = []
    vlines = []

    for start, end in intervals:
        s = start
        while s < end:
            row_start = (s // line_length) * line_length
            row_end = row_start + line_length
            fragment_end = min(end, row_end)
            y = base_y + (s // line_length) * ROW_SPACING

            # Rectangle
            rects.append(patches.Rectangle(
                (s - row_start, y - RECT_HEIGHT/2),
                fragment_end - s,
                RECT_HEIGHT
            ))

            # Vertical lines only at true atom boundaries
            if s == start:
                vlines.append([[s - row_start, y - VLINE_HEIGHT/2], [s - row_start, y + VLINE_HEIGHT/2]])
            if fragment_end == end:
                vlines.append([[fragment_end - row_start, y - VLINE_HEIGHT/2], [fragment_end - row_start, y + VLINE_HEIGHT/2]])

            s = fragment_end

    if rects:
        ax.add_collection(PatchCollection(rects, facecolor=color, edgecolor=color, alpha=1.0, zorder=3))
    if vlines:
        ax.add_collection(LineCollection(vlines, colors=color, linewidths=2, zorder=3))


def _draw_gap_lines(
    ax,
    intervals: list[tuple[int, int]],
    row: int,
    line_length: int,
    line_y: float,
) -> None:
    """
    Draw a genome line for a specific row, showing black segments only in the gaps
    where no atoms are present. This ensures that the line appears continuous except
    where predicted or true atoms exist in the given intervals.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The matplotlib axes to draw on.
    intervals : list of tuple[int, int]
        List of start-end positions of atoms for the genome.
    row : int
        The row index for the current line (used for wrapping by line_length).
    line_length : int
        Number of bases per row.
    line_y : float
        Y-coordinate on the plot where the line should be drawn.

    Returns
    -------
    None
    """
    row_start = row * line_length
    row_end = row_start + line_length
    intervals_in_row = []

    for start, end in intervals:
        if end <= row_start or start >= row_end:
            continue
        intervals_in_row.append((max(start, row_start), min(end, row_end)))

    intervals_in_row.sort()
    current = row_start
    lines = []

    for start, end in intervals_in_row:
        if start > current:
            lines.append([[current - row_start, line_y], [start - row_start, line_y]])
        current = end

    if current < row_end:
        lines.append([[current - row_start, line_y], [row_end - row_start, line_y]])

    if lines:
        ax.add_collection(LineCollection(lines, colors="black", linewidths=2, zorder=1))

# ---------------------------------------------------------------------
# Genome Atomization Visualization
# ---------------------------------------------------------------------

def plot_genome_atomization(
    genomes_file: Path,
    true_atoms_file: Path,
    predicted_atoms_file: Path,
    output_directory: Path,
) -> None:
    """
    Generate per-genome visualizations of predicted vs true atomization.

    Each genome gets its own PNG file showing the true atoms (blue) and predicted atoms (orange).

    Parameters
    ----------
    genomes_file : Path
        Input genomes FASTA file containing the genome sequences.
    true_atoms_file : Path
        Input GEESE TSV file with true (gold standard) atomization.
    predicted_atoms_file : Path
        Input GEESE TSV file with predicted atomization.
    output_directory : Path
        Directory where per-genome PNG visualizations are saved.

    Returns
    -------
    None

    Raises
    ------
    FileNotFoundError
        Raised if any of the input files do not exist.
    """
    if not genomes_file.is_file():
        raise FileNotFoundError(f"Genomes FASTA file not found: {genomes_file}")
    if not true_atoms_file.is_file():
        raise FileNotFoundError(f"True atoms file not found: {true_atoms_file}")
    if not predicted_atoms_file.is_file():
        raise FileNotFoundError(f"Predicted atoms file not found: {predicted_atoms_file}")

    output_directory.mkdir(parents=True, exist_ok=True)

    genome_dictionary = read_fasta(genomes_file)
    df_predicted = read_geese(predicted_atoms_file)
    df_true = read_geese(true_atoms_file)

    # Process each genome separately
    for genome_name, record in genome_dictionary.items():
        genome_length = len(record)
        line_length = max(10_000, genome_length // TARGET_ROWS)
        n_rows = (genome_length + line_length - 1) // line_length

        predicted_intervals = list(zip(
            df_predicted.loc[df_predicted["name"] == genome_name, "start"],
            df_predicted.loc[df_predicted["name"] == genome_name, "end"],
        ))

        true_intervals = list(zip(
            df_true.loc[df_true["name"] == genome_name, "start"],
            df_true.loc[df_true["name"] == genome_name, "end"],
        ))

        wrapped_true = _wrapped_intervals(true_intervals, line_length)
        wrapped_pred = _wrapped_intervals(predicted_intervals, line_length)

        fig, ax = plt.subplots(figsize=(12, 2 + 1.5 * n_rows))

        # Draw genome gap lines
        for row in range(n_rows):
            _draw_gap_lines(ax, true_intervals, row, line_length, 2 + row * ROW_SPACING)
            _draw_gap_lines(ax, predicted_intervals, row, line_length, 1 + row * ROW_SPACING)

        # Draw atoms
        _draw_atoms(ax, wrapped_true, line_length, 2, "tab:blue")
        _draw_atoms(ax, wrapped_pred, line_length, 1, "tab:orange")

        ax.set_xlim(0, line_length)
        ax.set_ylim(0.5, ROW_SPACING * n_rows)

        ticks = []
        labels = []

        for row in range(n_rows):
            ticks.append(1 + row * ROW_SPACING)
            labels.append("Predicted Atoms")
            ticks.append(2 + row * ROW_SPACING)
            labels.append("True Atoms")

        ax.set_yticks(ticks)
        ax.set_yticklabels(labels)

        ax.set_xlabel("Genome position")
        ax.set_title(f"Genome Atomization: {genome_name}")

        plt.tight_layout()
        output_file = output_directory / f"{genome_name}.png"
        plt.savefig(output_file)
        plt.close(fig)

        print(f"Saved visualization for genome '{genome_name}' ({genome_length} bp) -> {output_file}")
