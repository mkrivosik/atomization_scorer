"""
Tests for paf_to_geese() function.
"""

from pathlib import Path

import pytest

from atomization_scorer import paf_to_geese


# ------------------------------------------------------------------------------
# Helper: create minimal PAF file
# ------------------------------------------------------------------------------
def create_minimal_paf(tmp_path: Path):
    """Create a minimal PAF file for testing."""
    paf_file = tmp_path / "example.paf"
    paf_file.write_text(
        "query1\t1000\t0\t1000\t+\ttarget1|class_1\t1000\t0\t1000\t1000\t1000\t60\n"
        "query2\t2000\t0\t2000\t-\ttarget2|class_2\t2000\t100\t600\t500\t500\t50\n"
    )
    return paf_file


# ------------------------------------------------------------------------------
# Test: generates a valid GEESE file
# ------------------------------------------------------------------------------
def test_paf_to_geese_basic(tmp_path: Path):
    """paf_to_geese should convert a minimal PAF file to a valid GEESE TSV."""
    paf_file = create_minimal_paf(tmp_path)
    output_file = tmp_path / "example.geese"

    result = paf_to_geese(paf_file=paf_file, output_file=output_file)

    assert result == output_file
    assert output_file.is_file()

    lines = output_file.read_text().splitlines()
    assert len(lines) == 3
    assert lines[0] == "#name\tclass\tstart\tend"

    # Check first data line
    first_fields = lines[1].split("\t")
    assert first_fields[0] == "target1"
    assert first_fields[1] == "1"
    assert first_fields[2] == "0"
    assert first_fields[3] == "1000"

    # Check second data line
    second_fields = lines[2].split("\t")
    assert second_fields[0] == "target2"
    assert second_fields[1] == "2"
    assert second_fields[2] == "100"
    assert second_fields[3] == "600"


# ------------------------------------------------------------------------------
# Test: handles empty PAF file
# ------------------------------------------------------------------------------
def test_paf_to_geese_empty_file(tmp_path: Path):
    """paf_to_geese should create a GEESE file with only header if PAF is empty."""
    paf_file = tmp_path / "empty.paf"
    paf_file.write_text("")

    output_file = tmp_path / "empty.geese"

    result = paf_to_geese(paf_file=paf_file, output_file=output_file)

    assert result == output_file
    assert output_file.is_file()

    lines = output_file.read_text().splitlines()
    assert len(lines) == 1
    assert lines[0] == "#name\tclass\tstart\tend"


# ------------------------------------------------------------------------------
# Test: raises FileNotFoundError for missing PAF file
# ------------------------------------------------------------------------------
def test_paf_to_geese_missing_file(tmp_path: Path):
    """paf_to_geese should raise FileNotFoundError if input PAF file does not exist."""
    missing_file = tmp_path / "missing.paf"
    output_file = tmp_path / "output.geese"

    with pytest.raises(FileNotFoundError):
        paf_to_geese(paf_file=missing_file, output_file=output_file)
