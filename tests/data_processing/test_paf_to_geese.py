"""
Tests for the paf_to_geese () function.
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
# Helper: create custom PAF file
# ------------------------------------------------------------------------------
def create_custom_paf(tmp_path: Path, contents: str, filename: str = "example.paf"):
    """Create a custom PAF file for testing."""
    paf_file = tmp_path / filename
    paf_file.write_text(contents)
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

    # Check the first data line
    first_fields = lines[1].split("\t")
    assert first_fields[0] == "query1"
    assert first_fields[1] == "1"
    assert first_fields[2] == "0"
    assert first_fields[3] == "1000"

    # Check the second data line
    second_fields = lines[2].split("\t")
    assert second_fields[0] == "query2"
    assert second_fields[1] == "2"
    assert second_fields[2] == "0"
    assert second_fields[3] == "2000"


# ------------------------------------------------------------------------------
# Test: empty PAF file is converted to a header-only GEESE file
# ------------------------------------------------------------------------------
def test_paf_to_geese_empty_file(tmp_path: Path):
    """paf_to_geese should create a GEESE file with only a header if PAF is empty."""
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
# Test: missing PAF file is rejected
# ------------------------------------------------------------------------------
def test_paf_to_geese_missing_file(tmp_path: Path):
    """paf_to_geese should raise FileNotFoundError if the input PAF file does not exist."""
    missing_file = tmp_path / "missing.paf"
    output_file = tmp_path / "output.geese"

    with pytest.raises(FileNotFoundError):
        paf_to_geese(paf_file=missing_file, output_file=output_file)


# ------------------------------------------------------------------------------
# Test: malformed row with fewer than 6 fields is rejected
# ------------------------------------------------------------------------------
def test_paf_to_geese_malformed_row_raises_value_error(tmp_path: Path):
    """paf_to_geese should raise ValueError for a PAF row with fewer than 6 fields."""
    paf_file = create_custom_paf(
        tmp_path,
        "query1\t1000\t0\n",
        filename="malformed.paf",
    )
    output_file = tmp_path / "malformed.geese"

    with pytest.raises(ValueError, match="expected at least 6 tab-separated fields"):
        paf_to_geese(paf_file=paf_file, output_file=output_file)


# ------------------------------------------------------------------------------
# Test: missing class tag in the target header is rejected
# ------------------------------------------------------------------------------
def test_paf_to_geese_missing_class_tag_raises_value_error(tmp_path: Path):
    """paf_to_geese should raise ValueError if the target header lacks |class_."""
    paf_file = create_custom_paf(
        tmp_path,
        "query1\t1000\t0\t1000\t+\ttarget1\t1000\t0\t1000\t1000\t1000\t60\n",
        filename="missing_class.paf",
    )
    output_file = tmp_path / "missing_class.geese"

    with pytest.raises(ValueError, match=r"target header must contain '\|class_'"):
        paf_to_geese(paf_file=paf_file, output_file=output_file)


# ------------------------------------------------------------------------------
# Test: non-integer query start or end is rejected
# ------------------------------------------------------------------------------
def test_paf_to_geese_non_integer_start_or_end_raises_value_error(tmp_path: Path):
    """paf_to_geese should raise ValueError if the query start or end is not an integer."""
    invalid_rows = [
        (
            "query1\t1000\tstart\t1000\t+\ttarget1|class_1\t1000\t0\t1000\t1000\t1000\t60\n",
            "non_integer_start.paf",
            "expected an integer in field 3",
        ),
        (
            "query1\t1000\t0\tend\t+\ttarget1|class_1\t1000\t0\t1000\t1000\t1000\t60\n",
            "non_integer_end.paf",
            "expected an integer in field 4",
        ),
    ]

    for paf_contents, paf_name, error_match in invalid_rows:
        paf_file = create_custom_paf(tmp_path, paf_contents, filename=paf_name)
        output_file = tmp_path / f"{Path(paf_name).stem}.geese"

        with pytest.raises(ValueError, match=error_match):
            paf_to_geese(paf_file=paf_file, output_file=output_file)
