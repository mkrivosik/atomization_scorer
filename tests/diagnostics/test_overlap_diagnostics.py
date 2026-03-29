"""
Tests for the overlap_diagnostics.py module.
"""

import csv
import json
from pathlib import Path
from unittest.mock import patch

import pytest

from atomization_scorer import diagnose_paf_overlaps
from atomization_scorer.diagnostics import overlap_diagnostics as od


# --------------------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------------------
def _read_tsv_rows(path: Path) -> list[dict[str, str]]:
    """Read a TSV file into a list of row dictionaries."""
    header_aliases = {
        "query": "query_name",
        "anchor": "anchor_atom",
        "partner": "partner_atom",
        "a_q_start": "anchor_query_start",
        "a_q_end": "anchor_query_end",
        "p_q_start": "partner_query_start",
        "p_q_end": "partner_query_end",
        "a_t_start": "anchor_target_start",
        "a_t_end": "anchor_target_end",
        "p_t_start": "partner_target_start",
        "p_t_end": "partner_target_end",
        "a_strand": "anchor_strand",
        "p_strand": "partner_strand",
        "ov_start": "overlap_start",
        "ov_end": "overlap_end",
        "ov_len": "overlap_length",
        "ov_class": "overlap_class",
    }
    with path.open() as file:
        rows: list[dict[str, str]] = []
        for raw_row in csv.DictReader(file, delimiter="\t"):
            row: dict[str, str] = {}
            for key, value in raw_row.items():
                if key is None:
                    continue
                normalized_key = header_aliases.get(key.strip(), key.strip())
                row[normalized_key] = value.strip() if value is not None else ""
            if row and all(value and set(value) == {"-"} for value in row.values()):
                continue
            rows.append(row)
        return rows


def _read_json(path: Path):
    """Read a JSON file and return the parsed payload."""
    with path.open() as file:
        return json.load(file)


def _make_alignment(
    *,
    query_name: str = "genome1",
    query_start: int = 100,
    query_end: int = 400,
    strand: str = "+",
    target_name: str = "rep|A",
    target_start: int = 0,
    target_end: int = 300,
) -> od.DiagnosticAlignment:
    """Create a DiagnosticAlignment record for overlap-classification and detection tests."""
    return od.DiagnosticAlignment(
        query_name=query_name,
        query_length=5000,
        query_start=query_start,
        query_end=query_end,
        strand=strand,
        target_name=target_name,
        target_length=1000,
        target_start=target_start,
        target_end=target_end,
        matches=280,
        alignment_length=300,
        mapq=60,
        similarity=0.95,
    )


# --------------------------------------------------------------------------------------
# Test: parses valid PAF rows
# --------------------------------------------------------------------------------------
def test_parse_paf_alignment():
    """_parse_paf_alignment should parse valid mandatory and optional PAF fields."""
    alignment = od._parse_paf_alignment(
        "genome1\t5000\t100\t400\t-\trep|A\t600\t0\t300\t290\t300\t60\tdv:f:0.02\n",
        line_number=7,
    )

    assert alignment.query_name == "genome1"
    assert alignment.query_length == 5000
    assert alignment.query_start == 100
    assert alignment.query_end == 400
    assert alignment.strand == "-"
    assert alignment.target_name == "rep|A"
    assert alignment.target_length == 600
    assert alignment.target_start == 0
    assert alignment.target_end == 300
    assert alignment.matches == 290
    assert alignment.alignment_length == 300
    assert alignment.mapq == 60
    assert alignment.similarity == pytest.approx(0.98)


# --------------------------------------------------------------------------------------
# Test: malformed PAF rows are rejected
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("line", "expected_message"),
    [
        ("genome1\t5000\t100\n", "expected at least 12 tab-separated fields"),
        (
            "genome1\t5000\tstart\t400\t+\trep|A\t600\t0\t300\t290\t300\t60\n",
            "invalid numeric value in mandatory fields",
        ),
        (
            "genome1\t5000\t400\t400\t+\trep|A\t600\t0\t300\t290\t300\t60\n",
            "query start must be smaller than query end",
        ),
        (
            "genome1\t5000\t100\t400\t+\trep|A\t600\t300\t300\t290\t300\t60\n",
            "target start must be smaller than target end",
        ),
    ],
)
def test_parse_paf_alignment_rejects_malformed_rows(line: str, expected_message: str):
    """_parse_paf_alignment should raise ValueError for malformed PAF rows."""
    with pytest.raises(ValueError, match=expected_message):
        od._parse_paf_alignment(line, line_number=1)


# --------------------------------------------------------------------------------------
# Test: overlap classification labels edge-touching and internal cases
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("left", "right", "expected_class"),
    [
        (
            _make_alignment(query_start=100, query_end=400, target_name="rep|A"),
            _make_alignment(query_start=250, query_end=520, target_name="rep|B"),
            "both_edge",
        ),
        (
            _make_alignment(query_start=100, query_end=700, target_name="rep|A"),
            _make_alignment(query_start=250, query_end=400, target_name="rep|B"),
            "mixed_edge_internal",
        ),
        (
            _make_alignment(query_start=100, query_end=700, target_name="rep|A"),
            _make_alignment(query_start=200, query_end=600, target_name="rep|B"),
            "mixed_edge_internal",
        ),
    ],
)
def test_classify_overlap(
    left: od.DiagnosticAlignment,
    right: od.DiagnosticAlignment,
    expected_class: str,
):
    """_classify_overlap should return the expected overlap class for representative cases."""
    overlap_start = max(left.query_start, right.query_start)
    overlap_end = min(left.query_end, right.query_end)

    assert od._classify_overlap(overlap_start, overlap_end, left, right) == expected_class


# --------------------------------------------------------------------------------------
# Test: detects overlaps within one query and skips same-target pairs
# --------------------------------------------------------------------------------------
def test_detect_overlaps_returns_expected_records():
    """_detect_overlaps should return normalized overlap records within one query group."""
    overlaps = od._detect_overlaps(
        [
            _make_alignment(query_start=100, query_end=400, target_name="rep|B"),
            _make_alignment(query_start=250, query_end=520, target_name="rep|A"),
            _make_alignment(query_start=600, query_end=800, target_name="rep|C"),
            _make_alignment(query_start=260, query_end=380, target_name="rep|A"),
        ]
    )

    assert overlaps == [
        {
            "query_name": "genome1",
            "anchor_atom": "rep|A",
            "partner_atom": "rep|B",
            "anchor_query_start": "250",
            "anchor_query_end": "520",
            "partner_query_start": "100",
            "partner_query_end": "400",
            "anchor_atom_length": "1000",
            "partner_atom_length": "1000",
            "anchor_target_start": "0",
            "anchor_target_end": "300",
            "partner_target_start": "0",
            "partner_target_end": "300",
            "anchor_strand": "+",
            "partner_strand": "+",
            "overlap_start": "250",
            "overlap_end": "400",
            "overlap_length": "150",
            "overlap_class": "both_edge",
        },
        {
            "query_name": "genome1",
            "anchor_atom": "rep|A",
            "partner_atom": "rep|B",
            "anchor_query_start": "260",
            "anchor_query_end": "380",
            "partner_query_start": "100",
            "partner_query_end": "400",
            "anchor_atom_length": "1000",
            "partner_atom_length": "1000",
            "anchor_target_start": "0",
            "anchor_target_end": "300",
            "partner_target_start": "0",
            "partner_target_end": "300",
            "anchor_strand": "+",
            "partner_strand": "+",
            "overlap_start": "260",
            "overlap_end": "380",
            "overlap_length": "120",
            "overlap_class": "mixed_edge_internal",
        },
    ]


# --------------------------------------------------------------------------------------
# Test: missing input files are rejected
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("missing_side", "expected_message"),
    [
        ("paf", "PAF file not found"),
        ("representatives", "Representative FASTA file not found"),
    ],
)
def test_diagnose_paf_overlaps_missing_input_file(
    tmp_path: Path,
    missing_side: str,
    expected_message: str,
):
    """diagnose_paf_overlaps should raise FileNotFoundError if an input file is missing."""
    paf_file = tmp_path / "filtered.paf"
    representatives = tmp_path / "representatives.fa"
    paf_file.write_text("")
    representatives.write_text(">rep|A\nAAAA\n")

    if missing_side == "paf":
        paf_file = tmp_path / "missing.paf"
    else:
        representatives = tmp_path / "missing.fa"

    with pytest.raises(FileNotFoundError, match=expected_message):
        diagnose_paf_overlaps(
            paf_file=paf_file,
            representatives_fasta=representatives,
            output_directory=tmp_path / "overlap_diagnostics",
        )


# --------------------------------------------------------------------------------------
# Test: writes expected reports and anchor FASTA inputs
# --------------------------------------------------------------------------------------
@patch("atomization_scorer.diagnostics.dotter_runner.run_dotter_for_anchors")
def test_diagnose_paf_overlaps_writes_reports_and_anchor_fastas(
    mock_run_dotter,
    tmp_path: Path,
):
    """diagnose_paf_overlaps should emit overlap tables, summaries, and per-anchor FASTA files."""
    paf_file = tmp_path / "filtered.paf"
    paf_file.write_text(
        "genome1\t5000\t100\t400\t+\trep|A\t600\t0\t300\t290\t300\t60\tdv:f:0.02\n"
        "genome1\t5000\t250\t520\t-\trep|B\t700\t50\t320\t250\t270\t55\tdv:f:0.03\n"
        "genome1\t5000\t520\t700\t+\trep|C\t500\t0\t180\t180\t180\t60\tdv:f:0.00\n"
        "genome2\t5000\t0\t80\t+\trep|A\t600\t0\t80\t80\t80\t60\tdv:f:0.00\n"
        "genome2\t5000\t60\t140\t+\trep|D\t600\t0\t80\t80\t80\t60\tdv:f:0.00\n"
    )
    representatives = tmp_path / "representatives.fa"
    representatives.write_text(
        ">rep|A\nAACCGGTT\n"
        ">rep|B\nATGC\n"
        ">rep|C\nGATTACA\n"
        ">rep|D\nTTAA\n"
    )

    output_directory = tmp_path / "overlap_diagnostics"
    result = diagnose_paf_overlaps(
        paf_file=paf_file,
        representatives_fasta=representatives,
        output_directory=output_directory,
        minimum_report_overlap_length=50,
        minimum_plot_overlap_length=100,
    )

    assert result == output_directory
    mock_run_dotter.assert_called_once_with(
        anchors_directory=output_directory / "anchors",
        extra_args=None,
    )

    overlaps = _read_json(output_directory / "overlaps.json")
    assert overlaps == [
        {
            "query_name": "genome1",
            "anchor": {
                "atom": "rep|A",
                "representative_atom_length": 600,
                "query_start": 100,
                "query_end": 400,
                "target_start": 0,
                "target_end": 300,
                "strand": "+",
            },
            "partner": {
                "atom": "rep|B",
                "representative_atom_length": 700,
                "query_start": 250,
                "query_end": 520,
                "target_start": 50,
                "target_end": 320,
                "strand": "-",
            },
            "overlap": {
                "start": 250,
                "end": 400,
                "length": 150,
                "class": "both_edge",
            },
        }
    ]

    summary = {row["metric"]: row["value"] for row in _read_tsv_rows(output_directory / "summary.tsv")}
    assert summary == {
        "total_filtered_alignments": "5",
        "total_overlapping_pairs": "2",
        "pairs_le_50bp": "1",
        "pairs_gt_50bp": "1",
        "pairs_gt_100bp": "1",
        "pairs_both_edge": "2",
        "pairs_mixed_edge_internal": "0",
        "anchors_with_reported_overlaps": "1",
        "query_anchor_pairs_with_reported_overlaps": "1",
        "max_overlap_length": "150",
    }

    assert _read_tsv_rows(output_directory / "anchor_genome_summary.tsv") == [
        {
            "query_name": "genome1",
            "anchor_atom": "rep|A",
            "n_partners": "1",
            "n_overlaps": "1",
            "max_overlap_length": "150",
            "n_both_edge": "1",
            "n_mixed_edge_internal": "0",
        },
    ]

    anchor_a_directory = output_directory / "anchors" / "rep_A"
    assert (anchor_a_directory / "X.fasta").is_file()
    assert (anchor_a_directory / "Y.fasta").is_file()
    assert (anchor_a_directory / "pairs.tsv").is_file()
    assert "partner=rep|B" in (anchor_a_directory / "Y.fasta").read_text()
    assert "ATGC" in (anchor_a_directory / "Y.fasta").read_text()
    assert not (output_directory / "anchors" / "rep_B").exists()


# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
# Test: Dotter runner is skipped when run_dotter is disabled
# --------------------------------------------------------------------------------------
@patch("atomization_scorer.diagnostics.dotter_runner.run_dotter_for_anchors")
def test_diagnose_paf_overlaps_skips_dotter_when_disabled(
    mock_run_dotter,
    tmp_path: Path,
):
    """diagnose_paf_overlaps should not call the Dotter runner when run_dotter is False."""
    paf_file = tmp_path / "filtered.paf"
    paf_file.write_text(
        "genome1\t5000\t100\t400\t+\trep|A\t600\t0\t300\t290\t300\t60\tdv:f:0.02\n"
        "genome1\t5000\t250\t520\t+\trep|B\t700\t50\t320\t250\t270\t55\tdv:f:0.03\n"
    )
    representatives = tmp_path / "representatives.fa"
    representatives.write_text(
        ">rep|A\nAACCGGTT\n"
        ">rep|B\nATGC\n"
    )

    diagnose_paf_overlaps(
        paf_file=paf_file,
        representatives_fasta=representatives,
        output_directory=tmp_path / "overlap_diagnostics",
        minimum_report_overlap_length=50,
        minimum_plot_overlap_length=100,
        run_dotter=False,
    )

    mock_run_dotter.assert_not_called()


# --------------------------------------------------------------------------------------
# Test: optionally runs Dotter after generating anchor FASTA inputs
# --------------------------------------------------------------------------------------
@patch("atomization_scorer.diagnostics.dotter_runner.run_dotter_for_anchors")
def test_diagnose_paf_overlaps_optionally_runs_dotter(
    mock_run_dotter,
    tmp_path: Path,
):
    """diagnose_paf_overlaps should optionally run Dotter after writing anchor FASTA inputs."""
    paf_file = tmp_path / "filtered.paf"
    paf_file.write_text(
        "genome1\t5000\t100\t400\t+\trep|A\t600\t0\t300\t290\t300\t60\tdv:f:0.02\n"
        "genome1\t5000\t250\t520\t+\trep|B\t700\t50\t320\t250\t270\t55\tdv:f:0.03\n"
    )
    representatives = tmp_path / "representatives.fa"
    representatives.write_text(
        ">rep|A\nAACCGGTT\n"
        ">rep|B\nATGC\n"
    )

    output_directory = tmp_path / "overlap_diagnostics"
    diagnose_paf_overlaps(
        paf_file=paf_file,
        representatives_fasta=representatives,
        output_directory=output_directory,
        minimum_report_overlap_length=50,
        minimum_plot_overlap_length=100,
        run_dotter=True,
        dotter_extra_args=["-v"],
    )

    mock_run_dotter.assert_called_once_with(
        anchors_directory=output_directory / "anchors",
        extra_args=["-v"],
    )


# --------------------------------------------------------------------------------------
# Test: reverse mode mirrors anchor summaries and FASTA inputs
# --------------------------------------------------------------------------------------
@patch("atomization_scorer.diagnostics.dotter_runner.run_dotter_for_anchors")
def test_diagnose_paf_overlaps_optionally_includes_reverse_anchor_outputs(
    mock_run_dotter,
    tmp_path: Path,
):
    """diagnose_paf_overlaps should add reverse anchor outputs only when include_reverse is enabled."""
    paf_file = tmp_path / "filtered.paf"
    paf_file.write_text(
        "genome1\t5000\t100\t400\t+\trep|A\t600\t0\t300\t290\t300\t60\tdv:f:0.02\n"
        "genome1\t5000\t250\t520\t+\trep|B\t700\t50\t320\t250\t270\t55\tdv:f:0.03\n"
    )
    representatives = tmp_path / "representatives.fa"
    representatives.write_text(
        ">rep|A\nAACCGGTT\n"
        ">rep|B\nATGC\n"
    )

    output_directory = tmp_path / "overlap_diagnostics"
    diagnose_paf_overlaps(
        paf_file=paf_file,
        representatives_fasta=representatives,
        output_directory=output_directory,
        minimum_report_overlap_length=50,
        minimum_plot_overlap_length=100,
        include_reverse=True,
        run_dotter=True,
    )

    mock_run_dotter.assert_called_once_with(
        anchors_directory=output_directory / "anchors",
        extra_args=None,
    )
    assert (output_directory / "anchors" / "rep_A" / "X.fasta").is_file()
    assert (output_directory / "anchors" / "rep_B" / "X.fasta").is_file()
    assert _read_tsv_rows(output_directory / "anchor_genome_summary.tsv") == [
        {
            "query_name": "genome1",
            "anchor_atom": "rep|A",
            "n_partners": "1",
            "n_overlaps": "1",
            "max_overlap_length": "150",
            "n_both_edge": "1",
            "n_mixed_edge_internal": "0",
        },
        {
            "query_name": "genome1",
            "anchor_atom": "rep|B",
            "n_partners": "1",
            "n_overlaps": "1",
            "max_overlap_length": "150",
            "n_both_edge": "1",
            "n_mixed_edge_internal": "0",
        },
    ]


# --------------------------------------------------------------------------------------
# Test: empty or below-threshold overlaps still produce summary files
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("contents", "report_threshold", "plot_threshold", "expected_total_overlaps"),
    [
        (
            "genome1\t5000\t0\t100\t+\trep|A\t600\t0\t100\t100\t100\t60\tdv:f:0.00\n"
            "genome1\t5000\t100\t200\t+\trep|B\t600\t0\t100\t100\t100\t60\tdv:f:0.00\n",
            10,
            20,
            "0",
        ),
        (
            "genome1\t5000\t0\t80\t+\trep|A\t600\t0\t80\t80\t80\t60\tdv:f:0.00\n"
            "genome1\t5000\t60\t140\t+\trep|B\t600\t0\t80\t80\t80\t60\tdv:f:0.00\n",
            30,
            50,
            "1",
        ),
    ],
)
def test_diagnose_paf_overlaps_handles_empty_or_below_threshold_cases(
    tmp_path: Path,
    contents: str,
    report_threshold: int,
    plot_threshold: int,
    expected_total_overlaps: str,
):
    """diagnose_paf_overlaps should still write report files when no overlap is reportable or plottable."""
    paf_file = tmp_path / "filtered.paf"
    paf_file.write_text(contents)
    representatives = tmp_path / "representatives.fa"
    representatives.write_text(
        ">rep|A\nAAAA\n"
        ">rep|B\nCCCC\n"
    )

    output_directory = tmp_path / "overlap_diagnostics"
    diagnose_paf_overlaps(
        paf_file=paf_file,
        representatives_fasta=representatives,
        output_directory=output_directory,
        minimum_report_overlap_length=report_threshold,
        minimum_plot_overlap_length=plot_threshold,
    )

    summary = {row["metric"]: row["value"] for row in _read_tsv_rows(output_directory / "summary.tsv")}
    assert summary["total_overlapping_pairs"] == expected_total_overlaps

    assert _read_json(output_directory / "overlaps.json") == []
    assert _read_tsv_rows(output_directory / "anchor_genome_summary.tsv") == []
    assert list((output_directory / "anchors").iterdir()) == []
