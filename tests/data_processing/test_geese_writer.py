"""
Tests for the write_geese() function.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
from pathlib import Path

import pandas as pd
import pytest

from atomization_scorer import read_geese, write_geese


# --------------------------------------------------------------------------------------
# Helper: build a valid GEESE DataFrame
# --------------------------------------------------------------------------------------
def _build_geese_dataframe() -> pd.DataFrame:
    """Create a valid in-memory GEESE table for writer tests."""
    return pd.DataFrame(
        {
            "name": ["genome1", "genome2"],
            "class": [1, 2],
            "start": [0, 50],
            "end": [100, 150],
        }
    )


# --------------------------------------------------------------------------------------
# Test: basic GEESE writing
# --------------------------------------------------------------------------------------
def test_write_geese_basic(tmp_path: Path):
    """write_geese should write a valid GEESE tabular file."""
    output = tmp_path / "out.geese"
    df_atoms = _build_geese_dataframe()

    path = write_geese(df_atoms=df_atoms, output_path=output)

    assert path == output
    assert output.exists()

    lines = output.read_text().splitlines()
    assert lines[0] == "#name\tclass\tstart\tend"
    assert lines[1] == "genome1\t1\t0\t100"
    assert lines[2] == "genome2\t2\t50\t150"


# --------------------------------------------------------------------------------------
# Test: parent directory creation
# --------------------------------------------------------------------------------------
def test_write_geese_creates_parent_dir(output_dir: Path):
    """write_geese should create an output directory if it does not exist."""
    nested_directory = output_dir / "nested" / "folder"
    output = nested_directory / "atoms.geese"

    path = write_geese(df_atoms=_build_geese_dataframe(), output_path=output)

    assert output.exists()
    assert path == output


# --------------------------------------------------------------------------------------
# Test: writing an empty GEESE table
# --------------------------------------------------------------------------------------
def test_write_geese_empty_dataframe(tmp_path: Path):
    """write_geese should write headers when the GEESE table contains no rows."""
    output = tmp_path / "empty.geese"
    df_atoms = pd.DataFrame(columns=["name", "class", "start", "end"])

    path = write_geese(df_atoms=df_atoms, output_path=output)

    assert output.exists()
    assert path == output
    assert output.read_text() == "#name\tclass\tstart\tend\n"


# --------------------------------------------------------------------------------------
# Test: extra columns are preserved
# --------------------------------------------------------------------------------------
def test_write_geese_preserves_extra_columns(tmp_path: Path):
    """write_geese should preserve columns beyond the required GEESE fields."""
    output = tmp_path / "extra_columns.geese"
    df_atoms = pd.DataFrame(
        {
            "name": ["genome1"],
            "atom_nr": [1],
            "class": [3],
            "strand": ["+"],
            "start": [10],
            "end": [20],
            "score": [0.95],
            "note": ["first"],
        }
    )

    write_geese(df_atoms=df_atoms, output_path=output)
    result = read_geese(geese_file=output)

    assert list(result.columns) == ["name", "atom_nr", "class", "strand", "start", "end", "score", "note"]
    assert result.loc[0, "atom_nr"] == 1
    assert result.loc[0, "strand"] == "+"
    assert result.loc[0, "score"] == 0.95
    assert result.loc[0, "note"] == "first"


# --------------------------------------------------------------------------------------
# Test: write_geese output can be read back with read_geese
# --------------------------------------------------------------------------------------
def test_write_geese_round_trip_with_read_geese(tmp_path: Path):
    """write_geese should produce GEESE that round-trips correctly through read_geese."""
    output = tmp_path / "round_trip.geese"
    df_atoms = pd.DataFrame(
        {
            "name": ["genome1", "genome2"],
            "atom_nr": [1, 2],
            "class": [1, 2],
            "strand": ["+", "-"],
            "start": [0, 50],
            "end": [100, 150],
        }
    )

    write_geese(df_atoms=df_atoms, output_path=output)
    result = read_geese(geese_file=output)

    pd.testing.assert_frame_equal(result, df_atoms)


# --------------------------------------------------------------------------------------
# Test: specific missing required columns are reported
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("columns", "error_message"),
    [
        (["name", "start", "end"], "Missing required columns: class"),
        (["name", "class", "end"], "Missing required columns: start"),
        (["name", "class", "start"], "Missing required columns: end"),
        (["class", "start", "end"], "Missing required columns: name"),
        (["name", "class"], "Missing required columns: start, end"),
    ],
)
def test_write_geese_missing_columns(columns: list[str], error_message: str, tmp_path: Path):
    """write_geese should raise ValueError with exact missing required column names."""
    output = tmp_path / "broken.geese"
    df_atoms = pd.DataFrame({column: [] for column in columns})

    with pytest.raises(ValueError, match=error_message):
        write_geese(df_atoms=df_atoms, output_path=output)
