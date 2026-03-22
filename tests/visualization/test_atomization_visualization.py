"""
Tests for the plot_atomization() function.
"""

from pathlib import Path
from typing import cast

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PatchCollection
import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from PIL import Image

from atomization_scorer import plot_atomization
from atomization_scorer.visualization import atomization_visualization as av


# --------------------------------------------------------------------------------------
# Helper: create synthetic genome and atomization input files
# --------------------------------------------------------------------------------------
def _write_atomization_inputs(
    tmp_path: Path,
    genome_length: int,
    true_intervals: list[tuple[int, int]],
    predicted_intervals: list[tuple[int, int]],
    genome_name: str = "test_genome",
) -> tuple[Path, Path, Path]:
    """Create synthetic FASTA and GEESE files for atomization visualization tests."""
    sample_fasta = tmp_path / "sample.fa"
    true_geese = tmp_path / "true.geese"
    predicted_geese = tmp_path / "predicted.geese"

    sample_record = SeqRecord(Seq("A" * genome_length), id=genome_name)
    SeqIO.write([sample_record], sample_fasta, "fasta")

    df_true = pd.DataFrame(
        {
            "name": [genome_name] * len(true_intervals),
            "class": [f"T{index}" for index in range(len(true_intervals))],
            "start": [start for start, _ in true_intervals],
            "end": [end for _, end in true_intervals],
        }
    )
    df_true.to_csv(true_geese, sep="\t", index=False)

    df_predicted = pd.DataFrame(
        {
            "name": [genome_name] * len(predicted_intervals),
            "class": [f"P{index}" for index in range(len(predicted_intervals))],
            "start": [start for start, _ in predicted_intervals],
            "end": [end for _, end in predicted_intervals],
        }
    )
    df_predicted.to_csv(predicted_geese, sep="\t", index=False)

    return sample_fasta, true_geese, predicted_geese


# --------------------------------------------------------------------------------------
# Helper: collect line segments from one LineCollection
# --------------------------------------------------------------------------------------
def _segments_as_tuples(collection: LineCollection) -> list[tuple[tuple[float, float], tuple[float, float]]]:
    """Convert LineCollection segments to tuples for stable assertions."""
    return [
        ((float(segment[0][0]), float(segment[0][1])), (float(segment[1][0]), float(segment[1][1])))
        for segment in collection.get_segments()
    ]


# --------------------------------------------------------------------------------------
# Test: whole-genome atomization visualization generation
# --------------------------------------------------------------------------------------
def test_plot_atomization(
    output_dir: Path,
    tmp_path: Path,
):
    """plot_atomization should create one PNG file for the genome."""
    sample_fasta, true_geese, predicted_geese = _write_atomization_inputs(
        tmp_path=tmp_path,
        genome_length=20_000,
        true_intervals=[(0, 4900), (5400, 9999), (10000, 18000)],
        predicted_intervals=[(500, 5400), (5400, 10100), (10800, 17000)],
    )

    plot_atomization(
        genomes_file=sample_fasta,
        true_atoms_file=true_geese,
        predicted_atoms_file=predicted_geese,
        output_directory=output_dir,
        target_rows=4,
        min_bases_per_row=5_000,
        max_bases_per_row=5_000,
    )

    png_files = sorted(output_dir.glob("*.png"))
    assert len(png_files) == 1, "Expected exactly one atomization visualization PNG file."

    with Image.open(png_files[0]) as img:
        img.verify()


# --------------------------------------------------------------------------------------
# Test: _draw_baseline draws only uncovered genome segments across wrapped rows
# --------------------------------------------------------------------------------------
def test_draw_baseline_segments_across_rows():
    """_draw_baseline should draw one black line segment for each uncovered row-local gap."""
    fig, ax = plt.subplots()

    av._draw_baseline(
        ax=ax,
        intervals=[(0, 4900), (5400, 9999), (10000, 18000)],
        genome_length=20_000,
        line_length=5_000,
        track_y=1.0,
    )

    assert len(ax.collections) == 1
    baseline_collection = ax.collections[0]
    assert isinstance(baseline_collection, LineCollection)
    assert _segments_as_tuples(baseline_collection) == [
        ((4900.0, 1.0), (5000.0, 1.0)),
        ((0.0, 4.0), (400.0, 4.0)),
        ((4999.0, 4.0), (5000.0, 4.0)),
        ((3000.0, 10.0), (5000.0, 10.0)),
    ]

    plt.close(fig)


# --------------------------------------------------------------------------------------
# Test: _draw_baseline does not draw anything when the whole genome is covered
# --------------------------------------------------------------------------------------
def test_draw_baseline_no_segments_when_fully_covered():
    """_draw_baseline should not add a collection if there are no uncovered gaps."""
    fig, ax = plt.subplots()

    av._draw_baseline(
        ax=ax,
        intervals=[(0, 5000), (5000, 10000)],
        genome_length=10_000,
        line_length=5_000,
        track_y=1.0,
    )

    assert len(ax.collections) == 0
    plt.close(fig)


# --------------------------------------------------------------------------------------
# Test: _draw_atom_track draws only true atom boundaries, not wrap boundaries
# --------------------------------------------------------------------------------------
def test_draw_atom_track_true_boundaries_only():
    """_draw_atom_track should add boundary markers only at true atom starts and ends."""
    fig, ax = plt.subplots()

    av._draw_atom_track(
        ax=ax,
        intervals=[(3000, 13000)],
        line_length=5_000,
        track_y=1.0,
        color="tab:orange",
    )

    assert len(ax.collections) == 2
    assert isinstance(ax.collections[0], PatchCollection)
    assert isinstance(ax.collections[1], LineCollection)
    patch_collection = cast(PatchCollection, ax.collections[0])
    boundary_collection = cast(LineCollection, ax.collections[1])
    assert len(patch_collection.get_paths()) == 3
    assert _segments_as_tuples(boundary_collection) == [
        ((3000.0, 0.825), (3000.0, 1.175)),
        ((3000.0, 6.825), (3000.0, 7.175)),
    ]

    plt.close(fig)


# --------------------------------------------------------------------------------------
# Test: _draw_atom_track shows true boundaries at wrapped row edges
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("interval", "expected_segments"),
    [
        (
            (5000, 7000),
            [((0.0, 3.825), (0.0, 4.175)), ((2000.0, 3.825), (2000.0, 4.175))],
        ),
        (
            (2000, 5000),
            [((2000.0, 0.825), (2000.0, 1.175)), ((5000.0, 0.825), (5000.0, 1.175))],
        ),
    ],
)
def test_draw_atom_track_boundary_markers_at_row_edges(
    interval: tuple[int, int],
    expected_segments: list[tuple[tuple[float, float], tuple[float, float]]],
):
    """_draw_atom_track should preserve true boundaries that land exactly on row edges."""
    fig, ax = plt.subplots()

    av._draw_atom_track(
        ax=ax,
        intervals=[interval],
        line_length=5_000,
        track_y=1.0,
        color="tab:blue",
    )

    assert len(ax.collections) == 2
    boundary_collection = cast(LineCollection, ax.collections[1])
    assert _segments_as_tuples(boundary_collection) == expected_segments

    plt.close(fig)


# --------------------------------------------------------------------------------------
# Test: adaptive row width follows min/max per-row bounds
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("genome_length", "target_rows", "min_bases_per_row", "max_bases_per_row", "expected_xlim"),
    [
        (20_000, 20, 10_000, 250_000, 10_000.5),
        (10_000_000, 20, 10_000, 250_000, 250_000.5),
    ],
)
def test_plot_atomization_adaptive_row_width(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    genome_length: int,
    target_rows: int,
    min_bases_per_row: int,
    max_bases_per_row: int,
    expected_xlim: float,
):
    """plot_atomization should adapt wrapped row width from genome length and bounds."""
    sample_fasta, true_geese, predicted_geese = _write_atomization_inputs(
        tmp_path=tmp_path,
        genome_length=genome_length,
        true_intervals=[(0, min(5_000, genome_length))],
        predicted_intervals=[(500, min(5_500, genome_length))],
    )
    captured_axes = []

    original_subplots = av.plt.subplots

    def recording_subplots(*args, **kwargs):
        fig, ax = original_subplots(*args, **kwargs)
        captured_axes.append(ax)
        return fig, ax

    monkeypatch.setattr(av.plt, "subplots", recording_subplots)

    plot_atomization(
        genomes_file=sample_fasta,
        true_atoms_file=true_geese,
        predicted_atoms_file=predicted_geese,
        output_directory=tmp_path / "output",
        target_rows=target_rows,
        min_bases_per_row=min_bases_per_row,
        max_bases_per_row=max_bases_per_row,
    )

    assert len(captured_axes) == 1
    x_limits = captured_axes[0].get_xlim()
    assert x_limits[0] == -0.5
    assert x_limits[1] == expected_xlim
    assert captured_axes[0].get_yticklabels()[0].get_text().startswith("Predicted")


# --------------------------------------------------------------------------------------
# Test: one-row rendering for genomes shorter than min_bases_per_row
# --------------------------------------------------------------------------------------
def test_plot_atomization_short_genome_single_row(output_dir: Path, tmp_path: Path):
    """plot_atomization should render a short genome as a single wrapped row."""
    sample_fasta, true_geese, predicted_geese = _write_atomization_inputs(
        tmp_path=tmp_path,
        genome_length=3_000,
        true_intervals=[(0, 1500)],
        predicted_intervals=[(500, 2500)],
    )

    plot_atomization(
        genomes_file=sample_fasta,
        true_atoms_file=true_geese,
        predicted_atoms_file=predicted_geese,
        output_directory=output_dir,
    )

    png_files = sorted(output_dir.glob("*.png"))
    assert len(png_files) == 1
    with Image.open(png_files[0]) as img:
        img.verify()


# --------------------------------------------------------------------------------------
# Test: render succeeds when one track has no atoms
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("true_intervals", "predicted_intervals"),
    [
        ([], [(500, 5400)]),
        ([(0, 4900)], []),
        ([], []),
    ],
)
def test_plot_atomization_empty_tracks(
    output_dir: Path,
    tmp_path: Path,
    true_intervals: list[tuple[int, int]],
    predicted_intervals: list[tuple[int, int]],
):
    """plot_atomization should still generate a valid figure when one or both tracks are empty."""
    sample_fasta, true_geese, predicted_geese = _write_atomization_inputs(
        tmp_path=tmp_path,
        genome_length=20_000,
        true_intervals=true_intervals,
        predicted_intervals=predicted_intervals,
    )

    plot_atomization(
        genomes_file=sample_fasta,
        true_atoms_file=true_geese,
        predicted_atoms_file=predicted_geese,
        output_directory=output_dir,
        target_rows=4,
        min_bases_per_row=5_000,
        max_bases_per_row=5_000,
    )

    png_files = sorted(output_dir.glob("*.png"))
    assert len(png_files) == 1
    with Image.open(png_files[0]) as img:
        img.verify()


# --------------------------------------------------------------------------------------
# Test: missing input files raise FileNotFoundError
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("missing_name", "expected_message"),
    [
        ("genomes_file", "Genomes FASTA file not found"),
        ("true_atoms_file", "True atoms file not found"),
        ("predicted_atoms_file", "Predicted atoms file not found"),
    ],
)
def test_plot_atomization_missing_file_raises(
    tmp_path: Path,
    missing_name: str,
    expected_message: str,
):
    """plot_atomization should raise FileNotFoundError when an input file is missing."""
    sample_fasta, true_geese, predicted_geese = _write_atomization_inputs(
        tmp_path=tmp_path,
        genome_length=20_000,
        true_intervals=[(0, 4900)],
        predicted_intervals=[(500, 5400)],
    )
    file_map = {
        "genomes_file": sample_fasta,
        "true_atoms_file": true_geese,
        "predicted_atoms_file": predicted_geese,
    }
    file_map[missing_name].unlink()

    with pytest.raises(FileNotFoundError, match=expected_message):
        plot_atomization(
            genomes_file=sample_fasta,
            true_atoms_file=true_geese,
            predicted_atoms_file=predicted_geese,
            output_directory=tmp_path / "output",
        )


# --------------------------------------------------------------------------------------
# Test: invalid intervals raise ValueError
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("true_intervals", "predicted_intervals", "expected_message"),
    [
        ([(-1, 10)], [(0, 5)], "contains a negative coordinate"),
        ([(5, 5)], [(0, 5)], "must satisfy start < end"),
        ([(0, 25_000)], [(0, 5)], "ends outside genome length"),
        ([(0, 10), (5, 15)], [(0, 5)], "plotting will continue"),
    ],
)
def test_plot_atomization_invalid_or_overlapping_intervals(
    output_dir: Path,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
    true_intervals: list[tuple[int, int]],
    predicted_intervals: list[tuple[int, int]],
    expected_message: str,
):
    """plot_atomization should reject invalid intervals and log overlapping ones."""
    sample_fasta, true_geese, predicted_geese = _write_atomization_inputs(
        tmp_path=tmp_path,
        genome_length=20_000,
        true_intervals=true_intervals,
        predicted_intervals=predicted_intervals,
    )

    if expected_message == "plotting will continue":
        with caplog.at_level("WARNING"):
            plot_atomization(
                genomes_file=sample_fasta,
                true_atoms_file=true_geese,
                predicted_atoms_file=predicted_geese,
                output_directory=output_dir,
                target_rows=4,
                min_bases_per_row=5_000,
                max_bases_per_row=5_000,
            )
        assert "plotting will continue" in caplog.text
    else:
        with pytest.raises(ValueError, match=expected_message):
            plot_atomization(
                genomes_file=sample_fasta,
                true_atoms_file=true_geese,
                predicted_atoms_file=predicted_geese,
                output_directory=output_dir,
            )


# --------------------------------------------------------------------------------------
# Test: invalid adaptive configuration raises ValueError
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("target_rows", "min_bases_per_row", "max_bases_per_row", "expected_message"),
    [
        (0, 10_000, 250_000, "target_rows must be a positive integer"),
        (20, 0, 250_000, "min_bases_per_row must be a positive integer"),
        (20, 10_000, 0, "max_bases_per_row must be a positive integer"),
        (20, 20_000, 10_000, "min_bases_per_row must be less than or equal to max_bases_per_row"),
    ],
)
def test_plot_atomization_invalid_configuration_raises(
    tmp_path: Path,
    target_rows: int,
    min_bases_per_row: int,
    max_bases_per_row: int,
    expected_message: str,
):
    """plot_atomization should reject invalid adaptive wrapping configuration."""
    sample_fasta, true_geese, predicted_geese = _write_atomization_inputs(
        tmp_path=tmp_path,
        genome_length=20_000,
        true_intervals=[(0, 4900)],
        predicted_intervals=[(500, 5400)],
    )

    with pytest.raises(ValueError, match=expected_message):
        plot_atomization(
            genomes_file=sample_fasta,
            true_atoms_file=true_geese,
            predicted_atoms_file=predicted_geese,
            output_directory=tmp_path / "output",
            target_rows=target_rows,
            min_bases_per_row=min_bases_per_row,
            max_bases_per_row=max_bases_per_row,
        )
