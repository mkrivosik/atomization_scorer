"""
Tests for the dotplot_inputs.py module.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
import csv
from pathlib import Path

import pytest

from atomization_scorer.diagnostics import dotplot_inputs as di


# --------------------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------------------
def _read_tsv_rows(path: Path) -> list[dict[str, str]]:
    """Read a TSV file into a list of row dictionaries."""
    with path.open() as file:
        return list(csv.DictReader(file, delimiter="\t"))


def _example_overlap(
    *,
    query_name: str = "genome1",
    anchor_atom: str = "rep|A",
    partner_atom: str = "rep|B",
    overlap_length: str = "150",
) -> dict[str, str]:
    """Build one representative overlap record for dotplot-input tests."""
    return {
        "query_name": query_name,
        "anchor_atom": anchor_atom,
        "partner_atom": partner_atom,
        "anchor_query_start": "100",
        "anchor_query_end": "400",
        "partner_query_start": "250",
        "partner_query_end": "520",
        "anchor_target_start": "0",
        "anchor_target_end": "300",
        "partner_target_start": "50",
        "partner_target_end": "320",
        "anchor_strand": "+",
        "partner_strand": "-",
        "overlap_start": "250",
        "overlap_end": "400",
        "overlap_length": overlap_length,
        "overlap_class": "both_edge",
    }


# --------------------------------------------------------------------------------------
# Test: path sanitization keeps safe names and hashes empty results
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ("rep|A", "rep_A"),
        ("class.17-alpha", "class.17-alpha"),
        ("|||", "anchor_98c4b7d3"),
    ],
)
def test_sanitize_path_component(value: str, expected: str):
    """_sanitize_path_component should preserve safe names and hash empty sanitized values."""
    assert di._sanitize_path_component(value) == expected


# --------------------------------------------------------------------------------------
# Test: partner FASTA header stays concise and traceable
# --------------------------------------------------------------------------------------
def test_build_partner_header():
    """_build_partner_header should include only the partner atom identifier."""
    header = di._build_partner_header(_example_overlap())

    assert header == "partner=rep|B"


# --------------------------------------------------------------------------------------
# Test: default mode creates one anchor directory per forward overlap only
# --------------------------------------------------------------------------------------
def test_write_anchor_dotplot_fastas_builds_forward_only_anchor_directories_by_default(tmp_path: Path):
    """write_anchor_dotplot_fastas should only create the forward anchor directory by default."""
    representatives = tmp_path / "representatives.fa"
    representatives.write_text(
        ">rep|A\nAACCGGTT\n"
        ">rep|B\nATGC\n"
    )

    anchors_directory = tmp_path / "anchors"
    di.write_anchor_dotplot_fastas(
        overlaps=[_example_overlap()],
        representatives_fasta=representatives,
        output_directory=anchors_directory,
        minimum_overlap_length=100,
    )

    anchor_a = anchors_directory / "rep_A"
    anchor_b = anchors_directory / "rep_B"

    assert (anchor_a / "X.fasta").is_file()
    assert (anchor_a / "Y.fasta").is_file()
    assert (anchor_a / "pairs.tsv").is_file()
    assert not anchor_b.exists()

    assert "AACCGGTT" in (anchor_a / "Y.fasta").read_text()
    assert "ATGC" in (anchor_a / "X.fasta").read_text()
    assert ">partner=rep|B" in (anchor_a / "X.fasta").read_text()

    anchor_a_pairs = _read_tsv_rows(anchor_a / "pairs.tsv")
    assert anchor_a_pairs[0]["partner_atom"] == "rep|B"


# --------------------------------------------------------------------------------------
# Test: reverse mode creates both anchor directories with swapped pair records
# --------------------------------------------------------------------------------------
def test_write_anchor_dotplot_fastas_optionally_builds_bidirectional_anchor_directories(tmp_path: Path):
    """write_anchor_dotplot_fastas should create both anchor views when include_reverse is enabled."""
    representatives = tmp_path / "representatives.fa"
    representatives.write_text(
        ">rep|A\nAACCGGTT\n"
        ">rep|B\nATGC\n"
    )

    anchors_directory = tmp_path / "anchors"
    di.write_anchor_dotplot_fastas(
        overlaps=[_example_overlap()],
        representatives_fasta=representatives,
        output_directory=anchors_directory,
        minimum_overlap_length=100,
        include_reverse=True,
    )

    anchor_a = anchors_directory / "rep_A"
    anchor_b = anchors_directory / "rep_B"

    assert (anchor_a / "X.fasta").is_file()
    assert (anchor_a / "Y.fasta").is_file()
    assert (anchor_a / "pairs.tsv").is_file()
    assert (anchor_b / "X.fasta").is_file()
    assert (anchor_b / "Y.fasta").is_file()
    assert (anchor_b / "pairs.tsv").is_file()

    assert "AACCGGTT" in (anchor_b / "X.fasta").read_text()
    anchor_b_pairs = _read_tsv_rows(anchor_b / "pairs.tsv")
    assert anchor_b_pairs[0]["partner_atom"] == "rep|A"


# --------------------------------------------------------------------------------------
# Test: partner entries are sorted by overlap length, partner atom, and query name
# --------------------------------------------------------------------------------------
def test_write_anchor_dotplot_fastas_sorts_partner_entries_for_each_anchor(tmp_path: Path):
    """write_anchor_dotplot_fastas should sort partner entries deterministically for one anchor."""
    representatives = tmp_path / "representatives.fa"
    representatives.write_text(
        ">rep|A\nAAAA\n"
        ">rep|B\nCCCC\n"
        ">rep|C\nGGGG\n"
        ">rep|D\nTTTT\n"
    )

    overlaps = [
        _example_overlap(query_name="genome2", partner_atom="rep|C", overlap_length="120"),
        _example_overlap(query_name="genome1", partner_atom="rep|D", overlap_length="150"),
        _example_overlap(query_name="genome3", partner_atom="rep|B", overlap_length="150"),
    ]

    anchors_directory = tmp_path / "anchors"
    di.write_anchor_dotplot_fastas(
        overlaps=overlaps,
        representatives_fasta=representatives,
        output_directory=anchors_directory,
        minimum_overlap_length=100,
    )

    partner_headers = [
        line[1:]
        for line in (anchors_directory / "rep_A" / "X.fasta").read_text().splitlines()
        if line.startswith(">")
    ]
    assert partner_headers == [
        "partner=rep|B",
        "partner=rep|D",
        "partner=rep|C",
    ]


# --------------------------------------------------------------------------------------
# Test: reverse anchor records swap anchor and partner coordinates correctly
# --------------------------------------------------------------------------------------
def test_write_anchor_dotplot_fastas_swaps_coordinates_in_reverse_pairs_table(tmp_path: Path):
    """write_anchor_dotplot_fastas should swap query and target coordinates in the reverse anchor pair row."""
    representatives = tmp_path / "representatives.fa"
    representatives.write_text(
        ">rep|A\nAAAA\n"
        ">rep|B\nCCCC\n"
    )

    anchors_directory = tmp_path / "anchors"
    di.write_anchor_dotplot_fastas(
        overlaps=[_example_overlap()],
        representatives_fasta=representatives,
        output_directory=anchors_directory,
        minimum_overlap_length=100,
        include_reverse=True,
    )

    reverse_row = _read_tsv_rows(anchors_directory / "rep_B" / "pairs.tsv")[0]
    assert reverse_row["anchor_atom"] == "rep|B"
    assert reverse_row["partner_atom"] == "rep|A"
    assert reverse_row["anchor_query_start"] == "250"
    assert reverse_row["anchor_query_end"] == "520"
    assert reverse_row["partner_query_start"] == "100"
    assert reverse_row["partner_query_end"] == "400"
    assert reverse_row["anchor_target_start"] == "50"
    assert reverse_row["anchor_target_end"] == "320"
    assert reverse_row["partner_target_start"] == "0"
    assert reverse_row["partner_target_end"] == "300"


# --------------------------------------------------------------------------------------
# Test: threshold filtering excludes only overlaps strictly below the plotting threshold
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("overlap_length", "minimum_overlap_length", "expected_empty"),
    [
        ("100", 100, False),
        ("80", 100, True),
    ],
)
def test_write_anchor_dotplot_fastas_skips_non_plottable_overlaps(
    tmp_path: Path,
    overlap_length: str,
    minimum_overlap_length: int,
    expected_empty: bool,
):
    """write_anchor_dotplot_fastas should skip overlaps that are strictly below the plotting threshold."""
    representatives = tmp_path / "representatives.fa"
    representatives.write_text(
        ">rep|A\nAAAA\n"
        ">rep|B\nCCCC\n"
    )

    anchors_directory = tmp_path / "anchors"
    result = di.write_anchor_dotplot_fastas(
        overlaps=[_example_overlap(overlap_length=overlap_length)],
        representatives_fasta=representatives,
        output_directory=anchors_directory,
        minimum_overlap_length=minimum_overlap_length,
    )

    assert result == anchors_directory
    if expected_empty:
        assert list(anchors_directory.iterdir()) == []
    else:
        assert (anchors_directory / "rep_A" / "X.fasta").is_file()
        assert (anchors_directory / "rep_A" / "Y.fasta").is_file()


# --------------------------------------------------------------------------------------
# Test: hashed fallback directory name is used when sanitization removes all characters
# --------------------------------------------------------------------------------------
def test_write_anchor_dotplot_fastas_uses_hashed_fallback_directory_name(tmp_path: Path):
    """write_anchor_dotplot_fastas should use a deterministic hashed fallback for empty sanitized names."""
    representatives = tmp_path / "representatives.fa"
    representatives.write_text(
        ">|||\nAAAA\n"
        ">rep|B\nCCCC\n"
    )

    anchors_directory = tmp_path / "anchors"
    di.write_anchor_dotplot_fastas(
        overlaps=[_example_overlap(anchor_atom="|||")],
        representatives_fasta=representatives,
        output_directory=anchors_directory,
        minimum_overlap_length=100,
    )

    assert (anchors_directory / "anchor_98c4b7d3" / "X.fasta").is_file()


# --------------------------------------------------------------------------------------
# Test: missing representative sequences are rejected
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("overlap", "error_message"),
    [
        (_example_overlap(anchor_atom="rep|missing"), "Representative FASTA is missing anchor atom"),
        (_example_overlap(partner_atom="rep|missing"), "Representative FASTA is missing partner atom"),
    ],
)
def test_write_anchor_dotplot_fastas_missing_representative_raises(
    tmp_path: Path,
    overlap: dict[str, str],
    error_message: str,
):
    """write_anchor_dotplot_fastas should raise ValueError when an anchor or partner is absent from FASTA."""
    representatives = tmp_path / "representatives.fa"
    representatives.write_text(
        ">rep|A\nAAAA\n"
        ">rep|B\nCCCC\n"
    )

    with pytest.raises(ValueError, match=error_message):
        di.write_anchor_dotplot_fastas(
            overlaps=[overlap],
            representatives_fasta=representatives,
            output_directory=tmp_path / "anchors",
            minimum_overlap_length=100,
        )
