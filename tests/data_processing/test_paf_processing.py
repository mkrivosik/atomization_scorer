"""
Tests for filter_paf() function.
"""

from pathlib import Path

import pytest

from atomization_scorer import filter_paf


# ------------------------------------------------------------------------------
# Helper: create minimal PAF file
# ------------------------------------------------------------------------------
def create_minimal_paf(tmp_path: Path):
    """Create a minimal PAF file for testing."""
    paf_file = tmp_path / "example.paf"
    paf_file.write_text(
        "genome1\t1000\t0\t1000\t+\trepresentative1\t1000\t0\t1000\t1000\t1000\t60\tde:f:0.02\n"
        "genome1\t1000\t0\t1000\t+\trepresentative2\t1000\t0\t1000\t800\t800\t50\tde:f:0.10\n"
        "genome1\t1000\t0\t1000\t+\trepresentative3\t1000\t0\t1000\t400\t400\t40\tde:f:0.01\n"
    )
    return paf_file


# ------------------------------------------------------------------------------
# Test: filtering by similarity and alignment length
# ------------------------------------------------------------------------------
def test_filter_paf_basic(tmp_path: Path):
    """filter_paf should keep only alignments meeting similarity and length thresholds."""
    paf_file = create_minimal_paf(tmp_path)
    output_file = tmp_path / "filtered.paf"

    result = filter_paf(
        paf_file=paf_file,
        output_file=output_file,
        minimum_similarity=0.95,
        minimum_alignment_length=500
    )

    assert result == output_file
    assert output_file.is_file()

    lines = output_file.read_text().splitlines()
    assert len(lines) == 1
    assert lines[0].startswith("genome1") and "representative1" in lines[0]


# ------------------------------------------------------------------------------
# Test: no lines pass filters
# ------------------------------------------------------------------------------
def test_filter_paf_no_lines(tmp_path: Path):
    """filter_paf should return empty file if no alignments meet thresholds."""
    paf_file = create_minimal_paf(tmp_path)
    output_file = tmp_path / "filtered_none.paf"

    filter_paf(
        paf_file=paf_file,
        output_file=output_file,
        minimum_similarity=0.99,
        minimum_alignment_length=2000
    )

    lines = output_file.read_text().splitlines()
    assert len(lines) == 0


# ------------------------------------------------------------------------------
# Test: handles empty PAF file
# ------------------------------------------------------------------------------
def test_filter_paf_empty_file(tmp_path: Path):
    """filter_paf should create an empty filtered file if PAF is empty."""
    paf_file = tmp_path / "empty.paf"
    paf_file.write_text("")

    output_file = tmp_path / "filtered_empty.paf"

    result = filter_paf(
        paf_file=paf_file,
        output_file=output_file,
        minimum_similarity=0.95,
        minimum_alignment_length=500
    )

    assert result == output_file
    assert output_file.is_file()

    lines = output_file.read_text().splitlines()
    assert len(lines) == 0


# ------------------------------------------------------------------------------
# Test: raises FileNotFoundError for missing PAF file
# ------------------------------------------------------------------------------
def test_filter_paf_missing_file(tmp_path: Path):
    """filter_paf should raise FileNotFoundError if the PAF file does not exist."""
    missing_file = tmp_path / "missing.paf"
    output_file = tmp_path / "filtered.paf"

    with pytest.raises(FileNotFoundError):
        filter_paf(paf_file=missing_file, output_file=output_file)
