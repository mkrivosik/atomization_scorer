"""
atomization_visualization.py

Interactive atomization visualization of predicted versus true atomization.

Functions
---------
_build_atom_source_data   : Build glyph data for one true or predicted atom track.
_build_ribbon_source_data : Build connector polygons for matched true/predicted atoms.
_render_genome_html       : Render one standalone interactive HTML visualization.
plot_atomization          : Generate one interactive HTML visualization per genome.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
from __future__ import annotations
import json
import logging
import hashlib
from pathlib import Path
from typing import Any

from atomization_scorer import read_fasta, read_geese

from .plotting_utils import (
    AtomRecord,
    build_class_color_map,
    compute_initial_window,
    get_atoms_for_genome,
    normalize_output_format,
    pair_atoms,
    sanitize_output_stem,
)

# --------------------------------------------------------------------------------------
# Logging
# --------------------------------------------------------------------------------------
log = logging.getLogger(__name__)

# --------------------------------------------------------------------------------------
# Layout Constants
# --------------------------------------------------------------------------------------
TRUE_TRACK_Y = 1.1
PREDICTED_TRACK_Y = 0.0
ATOM_HALF_HEIGHT = 0.18
BASELINE_HALF_HEIGHT = 0.05
BASELINE_COLOR = "#6B7C93"
PAGE_STYLE = """
body {
    margin: 0;
    padding: 0;
    background: linear-gradient(180deg, #eef3f8 0%, #f7fafc 100%);
    color: #182433;
    font-family: "IBM Plex Sans", "Segoe UI", sans-serif;
}
.bk-root {
    margin: 28px auto;
    max-width: 1280px;
}
.app-shell {
    margin: 0 auto;
    max-width: 1360px;
    padding: 28px 20px 40px 20px;
}
.app-header {
    margin: 0 auto 18px auto;
    padding: 22px 24px;
    border: 1px solid #d8e2ec;
    border-radius: 18px;
    background: rgba(255, 255, 255, 0.94);
    box-shadow: 0 18px 50px rgba(30, 55, 90, 0.08);
}
.app-kicker {
    display: inline-block;
    padding: 6px 11px;
    border-radius: 999px;
    background: #dbeafe;
    color: #1d4f91;
    font-size: 12px;
    font-weight: 700;
    letter-spacing: 0.02em;
    text-transform: uppercase;
}
.app-header h1 {
    margin: 14px 0 6px 0;
    font-size: 30px;
    line-height: 1.15;
}
.app-meta {
    margin: 0;
    color: #43566b;
}
.app-note {
    margin: 12px 0 0 0;
    max-width: 980px;
    color: #32465d;
}
.bk-Column {
    margin: 0 auto;
}
.viz-box {
    margin: 0 auto;
    padding: 22px 24px 26px 24px;
    background: rgba(255, 255, 255, 0.94);
    border: 1px solid #d8e2ec;
    border-radius: 18px;
    box-shadow: 0 18px 50px rgba(30, 55, 90, 0.08);
}
.selection-panel {
    padding: 14px 16px;
    border: 1px solid #d8e2ec;
    border-radius: 14px;
    background: #f8fbfd;
    margin-top: 8px;
}
.selection-status-text {
    color: #43566b;
    font-size: 14px;
}
.bk-input-group {
    border-radius: 12px;
}
.selection-table-scroll {
    overflow-x: auto;
    border: 1px solid #d8e2ec;
    border-radius: 14px;
    background: #ffffff;
}
.selection-html-table {
    width: 100%;
    min-width: 1080px;
    border-collapse: collapse;
    font-size: 14px;
}
.selection-html-table th,
.selection-html-table td {
    padding: 10px 12px;
    border-bottom: 1px solid #e6edf4;
    text-align: left;
    vertical-align: top;
}
.selection-html-table thead th {
    background: #f5f9fc;
    color: #314558;
    font-size: 12px;
    font-weight: 700;
    letter-spacing: 0.01em;
    text-transform: uppercase;
}
.selection-html-table tbody tr {
    cursor: pointer;
}
.selection-html-table tbody tr:hover {
    background: #f8fbfd;
}
.selection-html-table tbody tr.is-selected {
    background: #e8f2fb;
}
.selection-html-table tbody tr:last-child td {
    border-bottom: none;
}
.selection-checkbox-cell {
    width: 44px;
    text-align: center;
}
.selection-checkbox {
    width: 16px;
    height: 16px;
    cursor: pointer;
}
.selection-table-empty {
    padding: 18px 16px;
    color: #52667b;
}
.selection-actions {
    display: flex;
    gap: 12px;
    margin-top: 14px;
    flex-wrap: wrap;
}
.selection-button {
    appearance: none;
    border: 1px solid #bfd0e0;
    border-radius: 12px;
    background: #f6f9fc;
    color: #213244;
    cursor: pointer;
    font: inherit;
    font-weight: 600;
    padding: 11px 16px;
    transition: background 120ms ease, border-color 120ms ease, transform 120ms ease;
}
.selection-button:hover {
    background: #edf4fa;
    border-color: #a8bfd4;
}
.selection-button:active {
    transform: translateY(1px);
}
"""


# --------------------------------------------------------------------------------------
# Internal Helpers
# --------------------------------------------------------------------------------------
def _atom_signature(atom: AtomRecord) -> tuple[str, int, int, int, str]:
    """
    Build a stable signature for one atom record.

    Parameters
    ----------
    atom : AtomRecord
        One normalized atom record.

    Returns
    -------
    tuple[str, int, int, int, str]
        Signature suitable for membership tests.
    """
    return (
        atom["class_id"],
        atom["atom_number"],
        atom["start"],
        atom["end"],
        atom["source"],
    )


def _build_connector_polygon(
    true_atom: AtomRecord,
    predicted_atom: AtomRecord,
) -> tuple[list[float], list[float]]:
    """
    Build a straight-sided connector polygon between a matched true/predicted atom pair.

    Parameters
    ----------
    true_atom : AtomRecord
        True atom record.
    predicted_atom : AtomRecord
        Predicted atom record.

    Returns
    -------
    tuple[list[float], list[float]]
        X and Y coordinates for one ribbon polygon.
    """
    true_y = TRUE_TRACK_Y - ATOM_HALF_HEIGHT
    predicted_y = PREDICTED_TRACK_Y + ATOM_HALF_HEIGHT
    polygon = [
        (float(true_atom["start"]), true_y),
        (float(true_atom["end"]), true_y),
        (float(predicted_atom["end"]), predicted_y),
        (float(predicted_atom["start"]), predicted_y),
    ]
    xs = [point[0] for point in polygon]
    ys = [point[1] for point in polygon]
    return xs, ys


def _build_gap_source_data(
    atoms: list[AtomRecord],
    genome_length: int,
    track_y: float,
    source_label: str,
) -> dict[str, list[Any]]:
    """
    Build thin baseline rectangles for uncovered genome segments on one track.

    Parameters
    ----------
    atoms : list[AtomRecord]
        Atom records for one source and one genome.
    genome_length : int
        Total genome length in bases.
    track_y : float
        Vertical center of the track.
    source_label : str
        Track label used in the data source.

    Returns
    -------
    dict[str, list[Any]]
        ColumnDataSource-compatible dictionary for uncovered genome gaps.
    """
    data = {
        "start": [],
        "end": [],
        "top": [],
        "bottom": [],
        "source": [],
        "length": [],
        "fill_color": [],
        "line_color": [],
    }
    def _append_gap_segment(gap_start: int, gap_end: int) -> None:
        data["start"].append(gap_start)
        data["end"].append(gap_end)
        data["top"].append(track_y + BASELINE_HALF_HEIGHT)
        data["bottom"].append(track_y - BASELINE_HALF_HEIGHT)
        data["source"].append(source_label)
        data["length"].append(gap_end - gap_start)
        data["fill_color"].append(BASELINE_COLOR)
        data["line_color"].append(BASELINE_COLOR)

    covered_end = 0

    for atom in sorted(atoms, key=lambda record: (record["start"], record["end"])):
        start = max(0, int(atom["start"]))
        end = min(genome_length, int(atom["end"]))
        if start > covered_end:
            _append_gap_segment(covered_end, start)
        covered_end = max(covered_end, end)

    if covered_end < genome_length:
        _append_gap_segment(covered_end, genome_length)

    return data


def _build_atom_source_data(
    atoms: list[AtomRecord],
    class_colors: dict[str, str],
    outline_color: str,
    status_by_signature: dict[tuple[str, int, int, int, str], str],
    track_y: float,
    genome_sequence: str,
) -> dict[str, list[Any]]:
    """
    Build glyph data for one true or predicted atom track.

    Parameters
    ----------
    atoms : list[AtomRecord]
        Atom records for one source and one genome.
    class_colors : dict[str, str]
        Mapping from class identifier to color.
    outline_color : str
        Outline color used for this track.
    status_by_signature : dict[tuple[str, int, int, int, str], str]
        Match status keyed by atom signature.
    track_y : float
        Vertical center of the track.

    Returns
    -------
    dict[str, list[Any]]
        ColumnDataSource-compatible dictionary.
    """
    data = {
        "start": [],
        "end": [],
        "top": [],
        "bottom": [],
        "source": [],
        "class_id": [],
        "atom_number": [],
        "atom_id": [],
        "length": [],
        "sequence": [],
        "fill_color": [],
        "line_color": [],
        "fill_alpha": [],
        "status": [],
    }

    for atom in atoms:
        signature = _atom_signature(atom)
        status = status_by_signature.get(signature, "unmatched")
        data["start"].append(atom["start"])
        data["end"].append(atom["end"])
        data["top"].append(track_y + ATOM_HALF_HEIGHT)
        data["bottom"].append(track_y - ATOM_HALF_HEIGHT)
        data["source"].append(atom["source"])
        data["class_id"].append(atom["class_id"])
        data["atom_number"].append(atom["atom_number"])
        data["atom_id"].append(atom["atom_id"])
        data["length"].append(atom["length"])
        data["sequence"].append(genome_sequence[atom["start"]:atom["end"]])
        data["fill_color"].append(class_colors[atom["class_id"]])
        data["line_color"].append(outline_color)
        data["fill_alpha"].append(0.95 if status == "matched" else 0.55)
        data["status"].append(status)

    return data


def _build_selection_row_source_data(
    matched_pairs: list[tuple[AtomRecord, AtomRecord]],
    unmatched_true: list[AtomRecord],
    unmatched_predicted: list[AtomRecord],
) -> dict[str, list[Any]]:
    """
    Build selection-table rows for matched and unmatched atoms.

    Parameters
    ----------
    matched_pairs : list[tuple[AtomRecord, AtomRecord]]
        Matched true/predicted atom pairs.
    unmatched_true : list[AtomRecord]
        True atoms without a matching predicted atom.
    unmatched_predicted : list[AtomRecord]
        Predicted atoms without a matching true atom.

    Returns
    -------
    dict[str, list[Any]]
        ColumnDataSource-compatible dictionary for selection-table rows.
    """
    data = {
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

    def _append_row(
        row_key: str,
        status: str,
        class_id: str,
        p_atom: AtomRecord | None,
        t_atom: AtomRecord | None,
    ) -> None:
        data["row_key"].append(row_key)
        data["status"].append(status)
        data["class_id"].append(class_id)
        data["predicted_atom_nr"].append(p_atom["atom_number"] if p_atom is not None else None)
        data["predicted_start"].append(p_atom["start"] if p_atom is not None else None)
        data["predicted_end"].append(p_atom["end"] if p_atom is not None else None)
        data["predicted_length"].append(p_atom["length"] if p_atom is not None else None)
        data["true_atom_nr"].append(t_atom["atom_number"] if t_atom is not None else None)
        data["true_start"].append(t_atom["start"] if t_atom is not None else None)
        data["true_end"].append(t_atom["end"] if t_atom is not None else None)
        data["true_length"].append(t_atom["length"] if t_atom is not None else None)

    for true_atom, predicted_atom in matched_pairs:
        _append_row(
            row_key=(
                f"matched|{true_atom['class_id']}|{true_atom['atom_number']}|"
                f"{true_atom['start']}|{true_atom['end']}|"
                f"{predicted_atom['atom_number']}|{predicted_atom['start']}|{predicted_atom['end']}"
            ),
            status="matched",
            class_id=true_atom["class_id"],
            p_atom=predicted_atom,
            t_atom=true_atom,
        )

    for true_atom in unmatched_true:
        _append_row(
            row_key=f"true_only|{true_atom['class_id']}|{true_atom['atom_number']}|{true_atom['start']}|{true_atom['end']}",
            status="missing predicted",
            class_id=true_atom["class_id"],
            p_atom=None,
            t_atom=true_atom,
        )

    for predicted_atom in unmatched_predicted:
        _append_row(
            row_key=(
                f"predicted_only|{predicted_atom['class_id']}|{predicted_atom['atom_number']}|"
                f"{predicted_atom['start']}|{predicted_atom['end']}"
            ),
            status="unexpected predicted",
            class_id=predicted_atom["class_id"],
            p_atom=predicted_atom,
            t_atom=None,
        )

    return data


def _build_empty_selection_table_data() -> dict[str, list[Any]]:
    """
    Build an empty selection-table data payload.

    Returns
    -------
    dict[str, list[Any]]
        Empty ColumnDataSource-compatible dictionary for the selected-atoms table.
    """
    return {
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


def _build_ribbon_source_data(
    matched_pairs: list[tuple[AtomRecord, AtomRecord]],
    class_colors: dict[str, str],
) -> dict[str, list[Any]]:
    """
    Build connector polygons for matched true/predicted atoms.

    Parameters
    ----------
    matched_pairs : list[tuple[AtomRecord, AtomRecord]]
        Matched true/predicted atom pairs.
    class_colors : dict[str, str]
        Mapping from class identifier to color.

    Returns
    -------
    dict[str, list[Any]]
        ColumnDataSource-compatible dictionary for ribbon polygons.
    """
    data = {
        "xs": [],
        "ys": [],
        "class_id": [],
        "atom_number": [],
        "true_start": [],
        "true_end": [],
        "predicted_start": [],
        "predicted_end": [],
        "fill_color": [],
        "line_color": [],
    }

    for true_atom, predicted_atom in matched_pairs:
        xs, ys = _build_connector_polygon(true_atom, predicted_atom)
        data["xs"].append(xs)
        data["ys"].append(ys)
        data["class_id"].append(true_atom["class_id"])
        data["atom_number"].append(true_atom["atom_number"])
        data["true_start"].append(true_atom["start"])
        data["true_end"].append(true_atom["end"])
        data["predicted_start"].append(predicted_atom["start"])
        data["predicted_end"].append(predicted_atom["end"])
        data["fill_color"].append(class_colors[true_atom["class_id"]])
        data["line_color"].append(class_colors[true_atom["class_id"]])

    return data


def _load_bokeh() -> dict[str, Any]:
    """
    Import Bokeh components lazily so helper tests can run without the dependency.

    Raises
    ------
    ModuleNotFoundError
        Raised if Bokeh is not installed in the runtime environment.

    Returns
    -------
    dict[str, Any]
        Mapping of Bokeh components used by the renderer.
    """
    try:
        from bokeh.embed import components
        from bokeh.events import SelectionGeometry, Tap
        from bokeh.layouts import column, row
        from bokeh.models import (
            BoxSelectTool,
            ColumnDataSource,
            CustomJS,
            Div,
            FixedTicker,
            HoverTool,
            NumeralTickFormatter,
            PanTool,
            Range1d,
            RangeSlider,
            ResetTool,
            SaveTool,
            TapTool,
            WheelZoomTool,
        )
        from bokeh.plotting import figure
        from bokeh.resources import INLINE
    except ModuleNotFoundError as error:
        raise ModuleNotFoundError(
            "Bokeh is required for interactive atomization visualization. "
            "Install the project with visualization dependencies enabled."
        ) from error

    return {
        "BoxSelectTool": BoxSelectTool,
        "INLINE": INLINE,
        "ColumnDataSource": ColumnDataSource,
        "CustomJS": CustomJS,
        "Div": Div,
        "FixedTicker": FixedTicker,
        "HoverTool": HoverTool,
        "NumeralTickFormatter": NumeralTickFormatter,
        "PanTool": PanTool,
        "Range1d": Range1d,
        "RangeSlider": RangeSlider,
        "ResetTool": ResetTool,
        "SaveTool": SaveTool,
        "SelectionGeometry": SelectionGeometry,
        "Tap": Tap,
        "TapTool": TapTool,
        "WheelZoomTool": WheelZoomTool,
        "column": column,
        "components": components,
        "figure": figure,
        "row": row,
    }


def _render_genome_html(
    genome_name: str,
    genome_length: int,
    genome_sequence: str,
    true_atoms: list[AtomRecord],
    predicted_atoms: list[AtomRecord],
    matched_pairs: list[tuple[AtomRecord, AtomRecord]],
    unmatched_true: list[AtomRecord],
    unmatched_predicted: list[AtomRecord],
    class_colors: dict[str, str],
    initial_window_bases: int,
    figure_width: float,
    true_color: str,
    predicted_color: str,
) -> str:
    """
    Render one standalone interactive HTML visualization.

    Parameters
    ----------
    genome_name : str
        Genome identifier.
    genome_length : int
        Total genome length in bases.
    genome_sequence : str
        Genome sequence used for atom sequence details.
    true_atoms : list[AtomRecord]
        True atoms for the genome.
    predicted_atoms : list[AtomRecord]
        Predicted atoms for the genome.
    matched_pairs : list[tuple[AtomRecord, AtomRecord]]
        Matched true/predicted atom pairs.
    unmatched_true : list[AtomRecord]
        True atoms with no predicted partner.
    unmatched_predicted : list[AtomRecord]
        Predicted atoms with no true partner.
    class_colors : dict[str, str]
        Mapping from class identifier to class color.
    initial_window_bases : int
        Initial visible genome window size.
    figure_width : float
        Requested figure width in inches, converted to an approximate viewport width.
    true_color : str
        Outline color for the true track.
    predicted_color : str
        Outline color for the predicted track.

    Returns
    -------
    str
        Standalone HTML document string.
    """
    bokeh = _load_bokeh()

    viewport_width = max(900, int(figure_width * 100))
    initial_window_end = min(genome_length, initial_window_bases)
    minimum_window = max(10, min(initial_window_bases, max(10, genome_length // 500)))

    matched_true_signatures = {
        _atom_signature(true_atom) for true_atom, _ in matched_pairs
    }
    matched_predicted_signatures = {
        _atom_signature(predicted_atom) for _, predicted_atom in matched_pairs
    }
    empty_selection_table_data = _build_empty_selection_table_data()
    empty_selection_table_data_json = json.dumps(empty_selection_table_data)
    default_selection_status_html = (
        "Use Box Select on the main plot to collect atoms touched by the selected genome interval."
    )
    default_atom_details_html = (
        "<div style='padding:14px 16px;border:1px solid #d8e2ec;border-radius:14px;background:#f8fbfd;'>"
        "<b>Atom details</b><br>"
        "Click a true or predicted atom to inspect its sequence and coordinates."
        "</div>"
    )
    remove_selection_help_html = (
        "Choose one or more rows in the table, then use Remove Selected Rows."
    )
    cleared_selection_summary_html = "Selection table cleared."
    empty_selection_table_html = (
        "<div class='selection-table-empty'>"
        "No atoms added yet. Use Box Select on the main plot to add atoms to the table."
        "</div>"
    )

    true_status = {
        _atom_signature(atom): ("matched" if _atom_signature(atom) in matched_true_signatures else "missing")
        for atom in true_atoms
    }
    predicted_status = {
        _atom_signature(atom): (
            "matched" if _atom_signature(atom) in matched_predicted_signatures else "unexpected"
        )
        for atom in predicted_atoms
    }

    true_source = bokeh["ColumnDataSource"](
        _build_atom_source_data(
            atoms=true_atoms,
            class_colors=class_colors,
            outline_color=true_color,
            status_by_signature=true_status,
            track_y=TRUE_TRACK_Y,
            genome_sequence=genome_sequence,
        )
    )
    true_gap_source = bokeh["ColumnDataSource"](
        _build_gap_source_data(
            atoms=true_atoms,
            genome_length=genome_length,
            track_y=TRUE_TRACK_Y,
            source_label="true_gap",
        )
    )
    predicted_source = bokeh["ColumnDataSource"](
        _build_atom_source_data(
            atoms=predicted_atoms,
            class_colors=class_colors,
            outline_color=predicted_color,
            status_by_signature=predicted_status,
            track_y=PREDICTED_TRACK_Y,
            genome_sequence=genome_sequence,
        )
    )
    predicted_gap_source = bokeh["ColumnDataSource"](
        _build_gap_source_data(
            atoms=predicted_atoms,
            genome_length=genome_length,
            track_y=PREDICTED_TRACK_Y,
            source_label="predicted_gap",
        )
    )
    ribbon_source = bokeh["ColumnDataSource"](
        _build_ribbon_source_data(matched_pairs=matched_pairs, class_colors=class_colors)
    )
    all_selection_rows_source = bokeh["ColumnDataSource"](
        _build_selection_row_source_data(
            matched_pairs=matched_pairs,
            unmatched_true=unmatched_true,
            unmatched_predicted=unmatched_predicted,
        )
    )
    selected_atoms_source = bokeh["ColumnDataSource"](empty_selection_table_data)

    pan_tool = bokeh["PanTool"](dimensions="width")
    wheel_zoom_tool = bokeh["WheelZoomTool"](dimensions="width")
    hover_tool = bokeh["HoverTool"]()
    box_select_tool = bokeh["BoxSelectTool"](dimensions="width")
    tap_tool = bokeh["TapTool"]()
    tap_tool.mode = "replace"
    reset_tool = bokeh["ResetTool"]()
    save_tool = bokeh["SaveTool"]()
    plot = bokeh["figure"](
        title=f"Atomization: {genome_name}",
        width=viewport_width,
        height=420,
        x_range=bokeh["Range1d"](0, initial_window_end),
        y_range=bokeh["Range1d"](-0.6, 1.7),
        tools=[pan_tool, wheel_zoom_tool, hover_tool, box_select_tool, tap_tool, reset_tool, save_tool],
        toolbar_location="above",
        toolbar_sticky=False,
    )
    plot.x_range.bounds = (0, genome_length)
    plot.x_range.min_interval = minimum_window
    plot.x_range.max_interval = genome_length
    plot.xaxis.axis_label = "Genome position"
    plot.xaxis.formatter = bokeh["NumeralTickFormatter"](format="0,0")
    plot.yaxis.ticker = bokeh["FixedTicker"](ticks=[PREDICTED_TRACK_Y, TRUE_TRACK_Y])
    plot.yaxis.major_label_overrides = {
        PREDICTED_TRACK_Y: "Predicted",
        TRUE_TRACK_Y: "True",
    }
    plot.ygrid.grid_line_color = None
    plot.xgrid.grid_line_alpha = 0.15
    plot.outline_line_color = "#444444"
    plot.border_fill_color = "#ffffff"
    plot.background_fill_color = "#fcfdff"
    plot.toolbar.logo = None
    plot.toolbar.active_drag = pan_tool
    plot.toolbar.active_scroll = None
    plot.toolbar.active_inspect = [hover_tool]
    plot.min_border_right = 28

    plot.quad(
        left="start",
        right="end",
        bottom="bottom",
        top="top",
        source=true_gap_source,
        fill_color="fill_color",
        fill_alpha=0.95,
        line_color="line_color",
        line_alpha=1.0,
        line_width=1,
    )
    plot.quad(
        left="start",
        right="end",
        bottom="bottom",
        top="top",
        source=predicted_gap_source,
        fill_color="fill_color",
        fill_alpha=0.95,
        line_color="line_color",
        line_alpha=1.0,
        line_width=1,
    )
    plot.patches(
        xs="xs",
        ys="ys",
        source=ribbon_source,
        fill_color="fill_color",
        fill_alpha=0.15,
        line_color="line_color",
        line_alpha=0.55,
        line_width=1.5,
        hatch_pattern="/",
        hatch_alpha=0.45,
        hatch_color="line_color",
    )
    true_renderer = plot.quad(
        left="start",
        right="end",
        bottom="bottom",
        top="top",
        source=true_source,
        fill_color="fill_color",
        fill_alpha="fill_alpha",
        line_color="line_color",
        line_width=2,
    )
    predicted_renderer = plot.quad(
        left="start",
        right="end",
        bottom="bottom",
        top="top",
        source=predicted_source,
        fill_color="fill_color",
        fill_alpha="fill_alpha",
        line_color="line_color",
        line_width=2,
    )
    hover_tool.renderers = [true_renderer, predicted_renderer]
    hover_tool.tooltips = [
        ("Track", "@source"),
        ("Class", "@class_id"),
        ("Atom number", "@atom_number{0,0}"),
        ("Start", "@start{0,0}"),
        ("End", "@end{0,0}"),
        ("Length", "@length{0,0}"),
        ("Status", "@status"),
    ]
    box_select_tool.renderers = [true_renderer, predicted_renderer]
    tap_tool.renderers = [true_renderer, predicted_renderer]

    overview_plot = bokeh["figure"](
        title="Genome overview",
        width=viewport_width,
        height=140,
        x_range=bokeh["Range1d"](0, genome_length),
        y_range=bokeh["Range1d"](-0.6, 1.7),
        tools="",
        toolbar_location=None,
    )
    overview_plot.x_range.bounds = (0, genome_length)
    overview_plot.yaxis.ticker = bokeh["FixedTicker"](ticks=[PREDICTED_TRACK_Y, TRUE_TRACK_Y])
    overview_plot.yaxis.major_label_overrides = {
        PREDICTED_TRACK_Y: "Predicted",
        TRUE_TRACK_Y: "True",
    }
    overview_plot.xaxis.axis_label = "Genome position"
    overview_plot.xaxis.formatter = bokeh["NumeralTickFormatter"](format="0,0")
    overview_plot.ygrid.grid_line_color = None
    overview_plot.xgrid.grid_line_alpha = 0.10
    overview_plot.outline_line_color = "#444444"
    overview_plot.background_fill_color = "#f6f9fc"
    overview_plot.border_fill_color = "#ffffff"
    overview_plot.quad(
        left="start",
        right="end",
        bottom="bottom",
        top="top",
        source=true_gap_source,
        fill_color="fill_color",
        fill_alpha=0.55,
        line_color=None,
    )
    overview_plot.quad(
        left="start",
        right="end",
        bottom="bottom",
        top="top",
        source=predicted_gap_source,
        fill_color="fill_color",
        fill_alpha=0.55,
        line_color=None,
    )
    overview_plot.quad(
        left="start",
        right="end",
        bottom="bottom",
        top="top",
        source=true_source,
        fill_color="fill_color",
        fill_alpha=0.35,
        line_color=None,
    )
    overview_plot.quad(
        left="start",
        right="end",
        bottom="bottom",
        top="top",
        source=predicted_source,
        fill_color="fill_color",
        fill_alpha=0.35,
        line_color=None,
    )

    window_slider = bokeh["RangeSlider"](
        title="Visible genome window",
        start=0,
        end=max(1, genome_length),
        value=(0, initial_window_end),
        step=max(1, genome_length // 1000),
        format="0,0",
        width=viewport_width,
    )
    window_slider.js_on_change(
        "value",
        bokeh["CustomJS"](
            args={"x_range": plot.x_range, "genome_length": genome_length},
            code="""
                const start = Math.max(0, cb_obj.value[0])
                const end = Math.min(genome_length, cb_obj.value[1])
                if (!(end > start)) {
                    return
                }
                x_range.start = start
                x_range.end = end
            """,
        ),
    )
    range_sync_callback = bokeh["CustomJS"](
        args={"slider": window_slider, "x_range": plot.x_range, "genome_length": genome_length},
        code="""
            const start = Math.max(0, Math.floor(x_range.start))
            const end = Math.min(genome_length, Math.ceil(x_range.end))
            if (!(end > start)) {
                return
            }
            if (x_range.start !== start) {
                x_range.start = start
            }
            if (x_range.end !== end) {
                x_range.end = end
            }
            slider.value = [start, end]
        """,
    )
    plot.x_range.js_on_change("start", range_sync_callback)
    plot.x_range.js_on_change("end", range_sync_callback)

    atom_details_div_id = f"atom-details-{selected_atoms_source.id}"
    selection_table_container_id = f"selected-atoms-table-{selected_atoms_source.id}"
    selection_status_div_id = f"selection-status-{selected_atoms_source.id}"
    remove_button_id = f"remove-selected-atoms-{selected_atoms_source.id}"
    clear_button_id = f"clear-selected-atoms-{selected_atoms_source.id}"
    selection_table_render_function = f"renderSelectedAtomsTable_{selected_atoms_source.id}"

    plot.js_on_event(
        bokeh["SelectionGeometry"],
        bokeh["CustomJS"](
            args={
                "all_rows": all_selection_rows_source,
                "source": selected_atoms_source,
                "true_source": true_source,
                "predicted_source": predicted_source,
            },
            code=f"""
                const geometry = cb_obj.geometry
                if (!geometry || geometry.x0 == null || geometry.x1 == null || cb_obj.final !== true) {{
                    return
                }}
                const raw_start = Math.floor(Math.min(geometry.x0, geometry.x1))
                const raw_end = Math.ceil(Math.max(geometry.x0, geometry.x1))
                const sel_start = Math.max(0, raw_start)
                const sel_end = Math.min({genome_length}, raw_end)
                if (!(sel_end > sel_start)) {{
                    return
                }}
                const nextData = {{}}
                for (const [field, values] of Object.entries(source.data)) {{
                    nextData[field] = [...values]
                }}
                const existingKeys = new Set(nextData.row_key)
                let added = 0
                let addedMinStart = null
                let addedMaxEnd = null
                for (let index = 0; index < all_rows.data.row_key.length; index += 1) {{
                    const key = all_rows.data.row_key[index]
                    if (existingKeys.has(key)) {{
                        continue
                    }}
                    const predictedStart = all_rows.data.predicted_start[index]
                    const predictedEnd = all_rows.data.predicted_end[index]
                    const trueStart = all_rows.data.true_start[index]
                    const trueEnd = all_rows.data.true_end[index]
                    const predictedIntersects = (
                        predictedStart != null && predictedEnd != null && predictedStart < sel_end && sel_start < predictedEnd
                    )
                    const trueIntersects = (
                        trueStart != null && trueEnd != null && trueStart < sel_end && sel_start < trueEnd
                    )
                    if (!(predictedIntersects || trueIntersects)) {{
                        continue
                    }}
                    for (const field of Object.keys(nextData)) {{
                        nextData[field].push(all_rows.data[field][index])
                    }}
                    existingKeys.add(key)
                    added += 1
                    if (predictedStart != null) {{
                        if (addedMinStart === null || predictedStart < addedMinStart) addedMinStart = predictedStart
                        if (addedMaxEnd === null || predictedEnd > addedMaxEnd) addedMaxEnd = predictedEnd
                    }}
                    if (trueStart != null) {{
                        if (addedMinStart === null || trueStart < addedMinStart) addedMinStart = trueStart
                        if (addedMaxEnd === null || trueEnd > addedMaxEnd) addedMaxEnd = trueEnd
                    }}
                }}
                source.data = nextData
                true_source.selected.indices = []
                predicted_source.selected.indices = []
                source.change.emit()
                const atomRange = (addedMinStart !== null && addedMaxEnd !== null)
                    ? " spanning [" + Number(addedMinStart).toLocaleString() + ";" + Number(addedMaxEnd).toLocaleString() + ")"
                    : ""
                const statusEl = document.getElementById('{selection_status_div_id}')
                if (statusEl) {{
                    statusEl.innerHTML = "Added " + added + " row(s)" + atomRange + ". Total rows: " + source.data.row_key.length + "."
                }}
            """,
        ),
    )
    tap_tool.callback = bokeh["CustomJS"](
        args={},
        code=f"""
            const source = cb_data.source
            if (!source) {{
                return
            }}
            const escapeHtml = (value) => String(value)
                .replaceAll('&', '&amp;')
                .replaceAll('<', '&lt;')
                .replaceAll('>', '&gt;')
            const indices = [...source.selected.indices]
            if (indices.length === 0) {{
                return
            }}
            const index = indices[indices.length - 1]
            const track = source.data.source[index]
            const classId = source.data.class_id[index]
            const atomNumber = source.data.atom_number[index]
            const start = source.data.start[index]
            const end = source.data.end[index]
            const length = source.data.length[index]
            const status = source.data.status[index]
            const sequence = escapeHtml(source.data.sequence[index])
            const el = document.getElementById({json.dumps(atom_details_div_id)})
            if (!el) {{
                return
            }}
            el.innerHTML = (
                "<div style='padding:14px 16px;border:1px solid #d8e2ec;border-radius:14px;background:#f8fbfd;'>"
                + "<b>Atom details</b><br>"
                + "<span><b>Track:</b> " + escapeHtml(track) + "</span><br>"
                + "<span><b>Class:</b> " + escapeHtml(classId) + "</span><br>"
                + "<span><b>Atom number:</b> " + Number(atomNumber).toLocaleString() + "</span><br>"
                + "<span><b>Start:</b> " + Number(start).toLocaleString() + "</span><br>"
                + "<span><b>End:</b> " + Number(end).toLocaleString() + "</span><br>"
                + "<span><b>Length:</b> " + Number(length).toLocaleString() + "</span><br>"
                + "<span><b>Status:</b> " + escapeHtml(status) + "</span><br>"
                + "<div style='margin-top:10px;'><b>Sequence</b></div>"
                + "<pre style='margin:8px 0 0 0;max-height:180px;overflow:auto;white-space:pre-wrap;word-break:break-all;"
                + "padding:10px 12px;border:1px solid #d8e2ec;border-radius:10px;background:#ffffff;'>"
                + sequence + "</pre></div>"
            )
        """,
    )
    plot.js_on_event(
        bokeh["Tap"],
        bokeh["CustomJS"](
            args={
                "true_source": true_source,
                "predicted_source": predicted_source,
            },
            code=f"""
                setTimeout(() => {{
                    if (true_source.selected.indices.length === 0 && predicted_source.selected.indices.length === 0) {{
                        const el = document.getElementById({json.dumps(atom_details_div_id)})
                        if (el) {{
                            el.innerHTML = {json.dumps(default_atom_details_html)}
                        }}
                    }}
                }}, 0)
            """,
        ),
    )
    selected_atoms_source.js_on_change(
        "data",
        bokeh["CustomJS"](
            code=f"window.{selection_table_render_function} && window.{selection_table_render_function}()",
        ),
    )

    layout = bokeh["column"](
        plot,
        window_slider,
        overview_plot,
    )
    script, div = bokeh["components"](layout)
    resources = bokeh["INLINE"]
    resources_html = resources.render()
    atom_details_html = f"""
    <div id="{atom_details_div_id}" style="margin-top:8px;">{default_atom_details_html}</div>
    """
    selection_table_controls_html = f"""
    <section class="selection-panel">
      <b>Selected Atoms</b><br>
      <span id="{selection_status_div_id}" class="selection-status-text">{default_selection_status_html}</span>
      <div style="margin-top:10px;">
        <div id="{selection_table_container_id}" class="selection-table-scroll">
          {empty_selection_table_html}
        </div>
        <div class="selection-actions">
          <button id="{remove_button_id}" class="selection-button" type="button">Remove Selected Rows</button>
          <button id="{clear_button_id}" class="selection-button" type="button">Clear Selected Atoms</button>
        </div>
      </div>
    </section>
    """
    selection_table_controls_script = f"""
  <script>
  (() => {{
    const tableContainerId = {json.dumps(selection_table_container_id)};
    const statusDivId = {json.dumps(selection_status_div_id)};
    const removeButtonId = {json.dumps(remove_button_id)};
    const clearButtonId = {json.dumps(clear_button_id)};
    const selectedSourceId = {json.dumps(selected_atoms_source.id)};
    const emptyData = {empty_selection_table_data_json};
    const emptyTableHtml = {json.dumps(empty_selection_table_html)};
    const removeSelectionHelpHtml = {json.dumps(remove_selection_help_html)};
    const clearedSelectionSummaryHtml = {json.dumps(cleared_selection_summary_html)};
    const selectedRowKeys = new Set();

    const escapeHtml = (value) => String(value)
      .replaceAll('&', '&amp;')
      .replaceAll('<', '&lt;')
      .replaceAll('>', '&gt;');

    const formatText = (value) => value == null ? '' : escapeHtml(value);
    const formatNumber = (value) => value == null ? '' : Number(value).toLocaleString();

    const cloneData = (data) => {{
      const nextData = {{}};
      for (const [field, values] of Object.entries(data)) {{
        nextData[field] = [...values];
      }}
      return nextData;
    }};

    const getDocument = () => window.Bokeh?.documents?.[0] ?? null;
    const getSelectedSource = () => getDocument()?.get_model_by_id(selectedSourceId) ?? null;

    const updateSummary = (html) => {{
      const el = document.getElementById(statusDivId);
      if (el) {{
        el.innerHTML = html;
      }}
    }};

    const syncSelectedRowKeys = (rowKeys) => {{
      const activeKeys = new Set(rowKeys);
      for (const rowKey of [...selectedRowKeys]) {{
        if (!activeKeys.has(rowKey)) {{
          selectedRowKeys.delete(rowKey);
        }}
      }}
    }};

    const renderTable = () => {{
      const container = document.getElementById(tableContainerId);
      const source = getSelectedSource();
      if (!container || !source) {{
        return;
      }}
      const rowKeys = source.data.row_key ?? [];
      syncSelectedRowKeys(rowKeys);
      if (rowKeys.length === 0) {{
        container.innerHTML = emptyTableHtml;
        return;
      }}
      const rows = rowKeys.map((rowKey, index) => {{
        const encodedKey = encodeURIComponent(rowKey);
        const checked = selectedRowKeys.has(rowKey);
        return (
          `<tr class="${{checked ? 'is-selected' : ''}}" data-row-key="${{encodedKey}}">`
          + `<td class="selection-checkbox-cell"><input class="selection-checkbox" type="checkbox" data-row-key="${{encodedKey}}" ${{checked ? 'checked' : ''}}></td>`
          + `<td>${{formatText(source.data.status[index])}}</td>`
          + `<td>${{formatText(source.data.class_id[index])}}</td>`
          + `<td>${{formatNumber(source.data.predicted_atom_nr[index])}}</td>`
          + `<td>${{formatNumber(source.data.predicted_start[index])}}</td>`
          + `<td>${{formatNumber(source.data.predicted_end[index])}}</td>`
          + `<td>${{formatNumber(source.data.predicted_length[index])}}</td>`
          + `<td>${{formatNumber(source.data.true_atom_nr[index])}}</td>`
          + `<td>${{formatNumber(source.data.true_start[index])}}</td>`
          + `<td>${{formatNumber(source.data.true_end[index])}}</td>`
          + `<td>${{formatNumber(source.data.true_length[index])}}</td>`
          + `</tr>`
        );
      }}).join('');
      container.innerHTML = (
        `<table class="selection-html-table">`
        + `<thead><tr>`
        + `<th class="selection-checkbox-cell">Pick</th>`
        + `<th>Status</th>`
        + `<th>Class</th>`
        + `<th>Pred Atom</th>`
        + `<th>Pred Start</th>`
        + `<th>Pred End</th>`
        + `<th>Pred Length</th>`
        + `<th>True Atom</th>`
        + `<th>True Start</th>`
        + `<th>True End</th>`
        + `<th>True Length</th>`
        + `</tr></thead>`
        + `<tbody>${{rows}}</tbody>`
        + `</table>`
      );
    }};

    const bindControls = () => {{
      const container = document.getElementById(tableContainerId);
      const removeButton = document.getElementById(removeButtonId);
      const clearButton = document.getElementById(clearButtonId);
      if (!container || !removeButton || !clearButton || container.dataset.bound === 'true') {{
        return;
      }}
      container.dataset.bound = 'true';

      container.addEventListener('change', (event) => {{
        const target = event.target;
        if (!(target instanceof HTMLInputElement) || !target.classList.contains('selection-checkbox')) {{
          return;
        }}
        const rowKey = decodeURIComponent(target.dataset.rowKey || '');
        if (!rowKey) {{
          return;
        }}
        if (target.checked) {{
          selectedRowKeys.add(rowKey);
        }} else {{
          selectedRowKeys.delete(rowKey);
        }}
        renderTable();
      }});

      container.addEventListener('click', (event) => {{
        const target = event.target;
        if (target instanceof HTMLInputElement) {{
          return;
        }}
        const row = target instanceof Element ? target.closest('tr[data-row-key]') : null;
        if (!row) {{
          return;
        }}
        const rowKey = decodeURIComponent(row.dataset.rowKey || '');
        if (!rowKey) {{
          return;
        }}
        if (selectedRowKeys.has(rowKey)) {{
          selectedRowKeys.delete(rowKey);
        }} else {{
          selectedRowKeys.add(rowKey);
        }}
        renderTable();
      }});

      removeButton.addEventListener('click', () => {{
        const source = getSelectedSource();
        if (!source) {{
          return;
        }}
        if (selectedRowKeys.size === 0) {{
          updateSummary(removeSelectionHelpHtml);
          return;
        }}
        const currentData = source.data;
        const nextData = {{}};
        for (const field of Object.keys(currentData)) {{
          nextData[field] = [];
        }}
        const previousCount = currentData.row_key.length;
        for (let index = 0; index < currentData.row_key.length; index += 1) {{
          const rowKey = currentData.row_key[index];
          if (selectedRowKeys.has(rowKey)) {{
            continue;
          }}
          for (const field of Object.keys(nextData)) {{
            nextData[field].push(currentData[field][index]);
          }}
        }}
        selectedRowKeys.clear();
        source.data = nextData;
        source.selected.indices = [];
        source.change.emit();
        updateSummary(
          "Removed " + (previousCount - nextData.row_key.length) + " row(s). "
          + nextData.row_key.length + " row(s) remain in the table."
        );
      }});

      clearButton.addEventListener('click', () => {{
        const source = getSelectedSource();
        if (!source) {{
          return;
        }}
        selectedRowKeys.clear();
        source.data = cloneData(emptyData);
        source.selected.indices = [];
        source.change.emit();
        updateSummary(clearedSelectionSummaryHtml);
      }});
    }};

    const initialize = () => {{
      if (!getSelectedSource() || !document.getElementById(tableContainerId)) {{
        window.setTimeout(initialize, 100);
        return;
      }}
      bindControls();
      renderTable();
    }};

    window.{selection_table_render_function} = renderTable;
    if (document.readyState === 'loading') {{
      document.addEventListener('DOMContentLoaded', initialize, {{ once: true }});
    }} else {{
      initialize();
    }}
  }})();
  </script>
    """
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Atomization: {genome_name}</title>
  {resources_html}
  <style>{PAGE_STYLE}</style>
  <style>
    .app-header, .viz-box {{ max-width: {viewport_width}px; }}
  </style>
</head>
<body>
  <main class="app-shell">
    <section class="app-header">
      <div class="app-kicker">Interactive Genome View</div>
      <h1>Genome: {genome_name}</h1>
      <p class="app-meta">
        Length: {genome_length} bp | Matching overlaps: {len(matched_pairs)} |
        Missing predicted atoms: {len(unmatched_true)} |
        Unexpected predicted atoms: {len(unmatched_predicted)}
      </p>
      <p class="app-note">
        Compare predicted and true atoms along the genome. Atoms of the same class that overlap are connected.
        Hover over a bar to inspect metadata. Drag horizontally to move across the genome, use the mouse wheel to
        zoom, or use the bottom range slider as a bounded left-right navigator.
      </p>
    </section>
    <div class="viz-box">
      {div}
      {atom_details_html}
      {selection_table_controls_html}
    </div>
  </main>
  {script}
  {selection_table_controls_script}
</body>
</html>"""


# --------------------------------------------------------------------------------------
# Atomization Visualization
# --------------------------------------------------------------------------------------
def plot_atomization(
    genomes_file: Path,
    true_atoms_file: Path,
    predicted_atoms_file: Path,
    output_directory: Path,
    figure_width: float = 12.0,
    target_rows: int = 20,
    min_bases_per_row: int = 10_000,
    max_bases_per_row: int = 250_000,
    true_color: str = "#2C7FB8",
    predicted_color: str = "#F28E2B",
    output_format: str = "html",
) -> None:
    """
    Generate one interactive HTML visualization per genome comparing predicted and true atoms.

    Parameters
    ----------
    genomes_file : Path
        Input genomes FASTA file containing the genome sequences.
    true_atoms_file : Path
        Input GEESE file containing the true atomization.
    predicted_atoms_file : Path
        Input GEESE file containing the predicted atomization.
    output_directory : Path
        Path to the output directory where HTML files are stored.
    figure_width : float, optional, default=12.0
        Approximate viewport width in inches, converted to pixels for HTML rendering.
    target_rows : int, optional, default=20
        Used to derive the initial visible genome window size.
    min_bases_per_row : int, optional, default=10_000
        Minimum initial visible window size in bases.
    max_bases_per_row : int, optional, default=250_000
        Maximum initial visible window size in bases.
    true_color : str, optional, default="#2C7FB8"
        Outline color used for the true atom track.
    predicted_color : str, optional, default="#F28E2B"
        Outline color used for the predicted atom track.
    output_format : str, optional, default="html"
        Interactive output format. Only HTML is supported.

    Raises
    ------
    FileNotFoundError
        Raised if genomes_file, true_atoms_file, or predicted_atoms_file do not exist.
    ValueError
        Raised if the initial window configuration is invalid or if the output format is unsupported.

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

    log.info(
        (
            "Generating interactive atomization visualizations from "
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
    used_output_stems = {}

    for genome_name, sequence in genome_dictionary.items():
        genome_length = len(sequence)
        initial_window_bases = compute_initial_window(
            genome_length=genome_length,
            target_rows=target_rows,
            min_bases_per_row=min_bases_per_row,
            max_bases_per_row=max_bases_per_row,
        )

        true_atoms = get_atoms_for_genome(
            df=df_true,
            genome_name=genome_name,
            genome_length=genome_length,
            label="True",
            source="true",
        )
        predicted_atoms = get_atoms_for_genome(
            df=df_predicted,
            genome_name=genome_name,
            genome_length=genome_length,
            label="Predicted",
            source="predicted",
        )
        matched_pairs, unmatched_true, unmatched_predicted = pair_atoms(
            true_atoms=true_atoms,
            predicted_atoms=predicted_atoms,
        )
        class_colors = build_class_color_map(
            [atom["class_id"] for atom in true_atoms + predicted_atoms]
        )

        log.info(
            (
                "Rendering interactive atomization visualization for genome=%s "
                "length=%s true_atoms=%s predicted_atoms=%s matched=%s "
                "missing=%s unexpected=%s initial_window=%s"
            ),
            genome_name,
            genome_length,
            len(true_atoms),
            len(predicted_atoms),
            len(matched_pairs),
            len(unmatched_true),
            len(unmatched_predicted),
            initial_window_bases,
        )

        html = _render_genome_html(
            genome_name=genome_name,
            genome_length=genome_length,
            genome_sequence=str(sequence),
            true_atoms=true_atoms,
            predicted_atoms=predicted_atoms,
            matched_pairs=matched_pairs,
            unmatched_true=unmatched_true,
            unmatched_predicted=unmatched_predicted,
            class_colors=class_colors,
            initial_window_bases=initial_window_bases,
            figure_width=figure_width,
            true_color=true_color,
            predicted_color=predicted_color,
        )
        output_stem = sanitize_output_stem(genome_name)
        previous_genome = used_output_stems.get(output_stem)
        if previous_genome is not None and previous_genome != genome_name:
            digest = hashlib.sha1(genome_name.encode("utf-8")).hexdigest()[:8]
            output_stem = f"{output_stem}_{digest}"
        used_output_stems[output_stem] = genome_name
        output_path = output_directory / output_stem
        output_path.with_suffix(f".{normalized_format}").write_text(html, encoding="utf-8")

        log.info(
            "Saved interactive atomization visualization for genome=%s to %s.%s",
            genome_name,
            output_path,
            normalized_format,
        )
