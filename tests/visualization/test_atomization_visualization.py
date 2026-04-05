"""
Tests for the interactive plot_atomization() function.
"""

from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import pytest

from atomization_scorer import plot_atomization
from atomization_scorer.visualization import atomization_visualization as av


# --------------------------------------------------------------------------------------
# Helper: create synthetic genome and atomization input files
# --------------------------------------------------------------------------------------
def _write_atomization_inputs(
    tmp_path: Path,
    genome_length: int,
    true_rows: list[dict[str, object]],
    predicted_rows: list[dict[str, object]],
    genome_name: str = "test_genome",
) -> tuple[Path, Path, Path]:
    """Create synthetic FASTA and GEESE files for interactive visualization tests."""
    sample_fasta = tmp_path / "sample.fa"
    true_geese = tmp_path / "true.geese"
    predicted_geese = tmp_path / "predicted.geese"

    sample_record = SeqRecord(Seq("A" * genome_length), id=genome_name)
    SeqIO.write([sample_record], sample_fasta, "fasta")

    pd.DataFrame(true_rows).to_csv(true_geese, sep="\t", index=False)
    pd.DataFrame(predicted_rows).to_csv(predicted_geese, sep="\t", index=False)
    return sample_fasta, true_geese, predicted_geese


# --------------------------------------------------------------------------------------
# Test: plot_atomization writes one HTML file per genome
# --------------------------------------------------------------------------------------
def test_plot_atomization_writes_html(output_dir: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    """plot_atomization should write one HTML file per genome using the interactive renderer."""
    sample_fasta, true_geese, predicted_geese = _write_atomization_inputs(
        tmp_path=tmp_path,
        genome_length=20_000,
        true_rows=[
            {"name": "test_genome", "class": "A", "start": 0, "end": 4900},
            {"name": "test_genome", "class": "B", "start": 5400, "end": 9999},
        ],
        predicted_rows=[
            {"name": "test_genome", "class": "A", "start": 500, "end": 5400},
            {"name": "test_genome", "class": "C", "start": 10_800, "end": 17_000},
        ],
    )
    captured: dict[str, object] = {}

    def fake_render_genome_html(**kwargs):
        captured.update(kwargs)
        return "<html><body>interactive-plot</body></html>"

    monkeypatch.setattr(av, "_render_genome_html", fake_render_genome_html)

    plot_atomization(
        genomes_file=sample_fasta,
        true_atoms_file=true_geese,
        predicted_atoms_file=predicted_geese,
        output_directory=output_dir,
    )

    html_files = sorted(output_dir.glob("*.html"))
    assert len(html_files) == 1
    assert html_files[0].read_text(encoding="utf-8") == "<html><body>interactive-plot</body></html>"
    assert captured["genome_name"] == "test_genome"
    assert captured["genome_length"] == 20_000
    assert len(captured["matched_pairs"]) == 1
    assert len(captured["unmatched_true"]) == 1
    assert len(captured["unmatched_predicted"]) == 1
    assert set(captured["class_colors"]) == {"A", "B", "C"}


# --------------------------------------------------------------------------------------
# Test: real HTML output embeds Bokeh resources and interactive controls
# --------------------------------------------------------------------------------------
def test_plot_atomization_real_html_contains_inline_resources_and_controls(output_dir: Path, tmp_path: Path):
    """plot_atomization should emit standalone HTML with inline Bokeh assets and controls."""
    sample_fasta, true_geese, predicted_geese = _write_atomization_inputs(
        tmp_path=tmp_path,
        genome_length=20_000,
        true_rows=[
            {"name": "test_genome", "class": "A", "start": 0, "end": 4900},
            {"name": "test_genome", "class": "B", "start": 5400, "end": 9999},
        ],
        predicted_rows=[
            {"name": "test_genome", "class": "A", "start": 500, "end": 5400},
            {"name": "test_genome", "class": "B", "start": 5400, "end": 10100},
        ],
    )

    plot_atomization(
        genomes_file=sample_fasta,
        true_atoms_file=true_geese,
        predicted_atoms_file=predicted_geese,
        output_directory=output_dir,
    )

    html_path = output_dir / "test_genome.html"
    html = html_path.read_text(encoding="utf-8")
    assert "Interactive Genome View" in html
    assert "Visible genome window" in html
    assert "Selected Atoms" in html
    assert "Atom details" in html
    assert "background: linear-gradient(180deg" in html
    assert "@start{0,0}" in html
    assert "@atom_number{0,0}" in html
    assert '"format":"0,0"' in html
    assert '"active_scroll":null' in html
    assert "cdn.bokeh.org" not in html


# --------------------------------------------------------------------------------------
# Test: initial visible window still respects target/min/max bounds
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("genome_length", "target_rows", "minimum", "maximum", "expected_window"),
    [
        (20_000, 20, 10_000, 250_000, 10_000),
        (10_000_000, 20, 10_000, 250_000, 250_000),
        (3_000, 20, 10_000, 250_000, 3_000),
    ],
)
def test_plot_atomization_forwards_initial_window_configuration(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    genome_length: int,
    target_rows: int,
    minimum: int,
    maximum: int,
    expected_window: int,
):
    """plot_atomization should pass the computed initial genome window to the renderer."""
    sample_fasta, true_geese, predicted_geese = _write_atomization_inputs(
        tmp_path=tmp_path,
        genome_length=genome_length,
        true_rows=[{"name": "test_genome", "class": "A", "start": 0, "end": min(500, genome_length)}],
        predicted_rows=[{"name": "test_genome", "class": "A", "start": 0, "end": min(500, genome_length)}],
    )
    calls: list[int] = []

    def fake_render_genome_html(**kwargs):
        calls.append(int(kwargs["initial_window_bases"]))
        return "<html></html>"

    monkeypatch.setattr(av, "_render_genome_html", fake_render_genome_html)

    plot_atomization(
        genomes_file=sample_fasta,
        true_atoms_file=true_geese,
        predicted_atoms_file=predicted_geese,
        output_directory=tmp_path / "output",
        target_rows=target_rows,
        min_bases_per_row=minimum,
        max_bases_per_row=maximum,
    )

    assert calls == [expected_window]


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
        true_rows=[{"name": "test_genome", "class": "A", "start": 0, "end": 4900}],
        predicted_rows=[{"name": "test_genome", "class": "A", "start": 500, "end": 5400}],
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
# Test: invalid configuration still raises before rendering
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("target_rows", "minimum", "maximum", "expected_message"),
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
    minimum: int,
    maximum: int,
    expected_message: str,
):
    """plot_atomization should reject invalid interactive window configuration."""
    sample_fasta, true_geese, predicted_geese = _write_atomization_inputs(
        tmp_path=tmp_path,
        genome_length=20_000,
        true_rows=[{"name": "test_genome", "class": "A", "start": 0, "end": 4900}],
        predicted_rows=[{"name": "test_genome", "class": "A", "start": 500, "end": 5400}],
    )

    with pytest.raises(ValueError, match=expected_message):
        plot_atomization(
            genomes_file=sample_fasta,
            true_atoms_file=true_geese,
            predicted_atoms_file=predicted_geese,
            output_directory=tmp_path / "output",
            target_rows=target_rows,
            min_bases_per_row=minimum,
            max_bases_per_row=maximum,
        )


# --------------------------------------------------------------------------------------
# Test: only HTML output is supported by the new visualization
# --------------------------------------------------------------------------------------
def test_plot_atomization_rejects_non_html_output_format(tmp_path: Path):
    """plot_atomization should reject static output formats from the previous renderer."""
    sample_fasta, true_geese, predicted_geese = _write_atomization_inputs(
        tmp_path=tmp_path,
        genome_length=20_000,
        true_rows=[{"name": "test_genome", "class": "A", "start": 0, "end": 4900}],
        predicted_rows=[{"name": "test_genome", "class": "A", "start": 500, "end": 5400}],
    )

    with pytest.raises(ValueError, match="Unsupported output format"):
        plot_atomization(
            genomes_file=sample_fasta,
            true_atoms_file=true_geese,
            predicted_atoms_file=predicted_geese,
            output_directory=tmp_path / "output",
            output_format="png",
        )


# --------------------------------------------------------------------------------------
# Test: Bokeh loader exposes the components used by the interactive renderer
# --------------------------------------------------------------------------------------
def test_load_bokeh_returns_required_components():
    """_load_bokeh should expose the Bokeh components needed by the renderer."""
    components = av._load_bokeh()

    assert "figure" in components
    assert "ColumnDataSource" in components
    assert "CustomJS" in components
    assert "components" in components
    assert "NumberFormatter" in components
    assert "NumeralTickFormatter" in components
    assert "TapTool" in components
    assert "Tap" in components


# --------------------------------------------------------------------------------------
# Test: matched connectors use straight polygon edges instead of curved Bezier paths
# --------------------------------------------------------------------------------------
def test_build_connector_polygon_uses_straight_edges():
    """_build_connector_polygon should connect matched atoms with one straight-sided quadrilateral."""
    true_atom = {
        "start": 100,
        "end": 220,
        "class_id": "A",
        "atom_number": 1,
        "source": "true",
    }
    predicted_atom = {
        "start": 140,
        "end": 280,
        "class_id": "A",
        "atom_number": 1,
        "source": "predicted",
    }

    xs, ys = av._build_connector_polygon(true_atom=true_atom, predicted_atom=predicted_atom)

    assert xs == [100.0, 220.0, 280.0, 140.0]
    assert ys == [
        av.TRUE_TRACK_Y - av.ATOM_HALF_HEIGHT,
        av.TRUE_TRACK_Y - av.ATOM_HALF_HEIGHT,
        av.PREDICTED_TRACK_Y + av.ATOM_HALF_HEIGHT,
        av.PREDICTED_TRACK_Y + av.ATOM_HALF_HEIGHT,
    ]


# --------------------------------------------------------------------------------------
# Test: uncovered genome gaps render as thin baseline rectangles between atoms
# --------------------------------------------------------------------------------------
def test_build_gap_source_data_returns_uncovered_segments():
    """_build_gap_source_data should emit merged uncovered intervals for one track."""
    atoms = [
        {
            "start": 10,
            "end": 20,
            "class_id": "A",
            "atom_number": 1,
            "source": "true",
        },
        {
            "start": 18,
            "end": 30,
            "class_id": "B",
            "atom_number": 1,
            "source": "true",
        },
        {
            "start": 40,
            "end": 45,
            "class_id": "C",
            "atom_number": 1,
            "source": "true",
        },
    ]

    gap_data = av._build_gap_source_data(
        atoms=atoms,
        genome_length=50,
        track_y=av.TRUE_TRACK_Y,
        source_label="true_gap",
    )

    assert gap_data["start"] == [0, 30, 45]
    assert gap_data["end"] == [10, 40, 50]
    assert gap_data["source"] == ["true_gap", "true_gap", "true_gap"]
    assert gap_data["length"] == [10, 10, 5]
    assert gap_data["top"] == [av.TRUE_TRACK_Y + av.BASELINE_HALF_HEIGHT] * 3
    assert gap_data["bottom"] == [av.TRUE_TRACK_Y - av.BASELINE_HALF_HEIGHT] * 3
    assert gap_data["fill_color"] == [av.BASELINE_COLOR, av.BASELINE_COLOR, av.BASELINE_COLOR]
    assert gap_data["line_color"] == [av.BASELINE_COLOR, av.BASELINE_COLOR, av.BASELINE_COLOR]


# --------------------------------------------------------------------------------------
# Test: atom glyph data includes the sequence slice for click-driven details
# --------------------------------------------------------------------------------------
def test_build_atom_source_data_includes_atom_sequence():
    """_build_atom_source_data should carry the genome subsequence for each atom."""
    atoms = [
        {
            "start": 2,
            "end": 6,
            "class_id": "A",
            "atom_number": 9,
            "atom_id": "A:9",
            "length": 4,
            "source": "true",
        },
    ]

    atom_data = av._build_atom_source_data(
        atoms=atoms,
        class_colors={"A": "#123456"},
        outline_color="#000000",
        status_by_signature={("A", 9, 2, 6, "true"): "matched"},
        track_y=av.TRUE_TRACK_Y,
        genome_sequence="AACCGGTT",
    )

    assert atom_data["sequence"] == ["CCGG"]


# --------------------------------------------------------------------------------------
# Test: selection-row data captures matched and unmatched atom comparisons
# --------------------------------------------------------------------------------------
def test_build_selection_row_source_data_captures_comparison_rows():
    """_build_selection_row_source_data should emit rows for matched and unmatched atom comparisons."""
    matched_true = {
        "class_id": "A",
        "atom_number": 1,
        "start": 10,
        "end": 20,
        "length": 10,
        "source": "true",
    }
    matched_predicted = {
        "class_id": "A",
        "atom_number": 7,
        "start": 12,
        "end": 22,
        "length": 10,
        "source": "predicted",
    }
    unmatched_true = {
        "class_id": "B",
        "atom_number": 2,
        "start": 30,
        "end": 40,
        "length": 10,
        "source": "true",
    }
    unmatched_predicted = {
        "class_id": "C",
        "atom_number": 8,
        "start": 50,
        "end": 70,
        "length": 20,
        "source": "predicted",
    }

    row_data = av._build_selection_row_source_data(
        matched_pairs=[(matched_true, matched_predicted)],
        unmatched_true=[unmatched_true],
        unmatched_predicted=[unmatched_predicted],
    )

    assert row_data["status"] == ["matched", "missing predicted", "unexpected predicted"]
    assert row_data["class_id"] == ["A", "B", "C"]
    assert row_data["predicted_atom_nr"] == [7, None, 8]
    assert row_data["true_atom_nr"] == [1, 2, None]


# --------------------------------------------------------------------------------------
# Test: empty selection-table payload stays aligned with the rendered DataTable schema
# --------------------------------------------------------------------------------------
def test_build_empty_selection_table_data_returns_all_expected_columns():
    """_build_empty_selection_table_data should expose every selected-atoms table column."""
    assert av._build_empty_selection_table_data() == {
        "row_key": [],
        "status": [],
        "class_id": [],
        "predicted_atom_nr": [],
        "predicted_start": [],
        "predicted_end": [],
        "predicted_length": [],
        "true_atom_nr": [],
        "true_start": [],
        "true_end": [],
        "true_length": [],
    }


# --------------------------------------------------------------------------------------
# Test: rendered HTML keeps the interaction callback guards for tap/table behavior
# --------------------------------------------------------------------------------------
def test_plot_atomization_real_html_contains_interaction_regression_guards(
    output_dir: Path,
    tmp_path: Path,
):
    """plot_atomization should embed the JS guards used for repeated tap and table interactions."""
    sample_fasta, true_geese, predicted_geese = _write_atomization_inputs(
        tmp_path=tmp_path,
        genome_length=20_000,
        true_rows=[
            {"name": "test_genome", "class": "A", "start": 0, "end": 4900},
            {"name": "test_genome", "class": "B", "start": 5400, "end": 9999},
        ],
        predicted_rows=[
            {"name": "test_genome", "class": "A", "start": 500, "end": 5400},
            {"name": "test_genome", "class": "B", "start": 5400, "end": 10100},
        ],
    )

    plot_atomization(
        genomes_file=sample_fasta,
        true_atoms_file=true_geese,
        predicted_atoms_file=predicted_geese,
        output_directory=output_dir,
    )

    html = (output_dir / "test_genome.html").read_text(encoding="utf-8")
    assert "geometry.x0 == null || geometry.x1 == null || cb_obj.final !== true" in html
    assert "const indices = [...source.selected.indices]" in html
    assert "source.data = nextData" in html
    assert "true_source.selected.indices = []" in html
    assert '"mode":"replace"' in html
    assert "selectedRowKeys.clear()" in html
    assert "renderSelectedAtomsTable_" in html


# --------------------------------------------------------------------------------------
# Test: sanitized filename collisions still produce one HTML per genome
# --------------------------------------------------------------------------------------
def test_plot_atomization_avoids_output_filename_collisions(tmp_path: Path):
    """plot_atomization should not let sanitized genome names overwrite each other."""
    sample_fasta, true_geese, predicted_geese = _write_atomization_inputs(
        tmp_path=tmp_path,
        genome_length=500,
        true_rows=[
            {"name": "genome/1", "class": "A", "start": 0, "end": 100},
            {"name": "genome:1", "class": "B", "start": 20, "end": 120},
        ],
        predicted_rows=[
            {"name": "genome/1", "class": "A", "start": 5, "end": 105},
            {"name": "genome:1", "class": "B", "start": 25, "end": 125},
        ],
        genome_name="genome/1",
    )
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    SeqIO.write(
        [
            SeqRecord(Seq("A" * 500), id="genome/1"),
            SeqRecord(Seq("C" * 500), id="genome:1"),
        ],
        sample_fasta,
        "fasta",
    )

    plot_atomization(
        genomes_file=sample_fasta,
        true_atoms_file=true_geese,
        predicted_atoms_file=predicted_geese,
        output_directory=tmp_path / "output",
    )

    html_files = sorted((tmp_path / "output").glob("*.html"))
    assert len(html_files) == 2
