"""
Tests for the read_geese() function.
"""

from pathlib import Path

import pandas as pd
import pytest

from atomization_scorer import read_geese


# ---------------------------------------------------------------------
# Test: loads a valid GEESE tabular file and contains all required columns
# ---------------------------------------------------------------------
def test_read_geese_valid(tmp_path: Path):
    """read_geese should load a valid GEESE tabular file and contain all required columns."""
    geese_file = tmp_path / "example.geese"
    geese_file.write_text(
        "#name\tatom_nr\tclass\tstrand\tstart\tend\n"
        "gene1\t1\t1\t+\t0\t100\n"
        "gene2\t2\t2\t-\t50\t150\n"
    )

    df_atoms = read_geese(geese_file=geese_file)

    required_columns = ['name', 'atom_nr', 'class', 'strand', 'start', 'end']
    for column in required_columns:
        assert column in df_atoms.columns

    assert df_atoms.shape[0] == 2
    assert df_atoms.loc[0, 'name'] == "gene1"
    assert df_atoms.loc[1, 'atom_nr'] == 2


# ---------------------------------------------------------------------
# Test: '#name' column is renamed to 'name'
# ---------------------------------------------------------------------
def test_read_geese_renames_hash_name_column(tmp_path: Path):
    """read_geese should rename a '#name' column to 'name'."""
    geese_file = tmp_path / "rename.geese"
    geese_file.write_text(
        "#name\tclass\tstart\tend\n"
        "gene1\t1\t0\t100\n"
    )

    df_atoms = read_geese(geese_file=geese_file)

    assert "name" in df_atoms.columns
    assert "#name" not in df_atoms.columns
    assert df_atoms.loc[0, "name"] == "gene1"


# ---------------------------------------------------------------------
# Test: extra columns are preserved
# ---------------------------------------------------------------------
def test_read_geese_preserves_extra_columns(tmp_path: Path):
    """read_geese should preserve columns beyond the required GEESE fields."""
    geese_file = tmp_path / "extra_columns.geese"
    geese_file.write_text(
        "#name\tatom_nr\tclass\tstrand\tstart\tend\tscore\tnote\n"
        "gene1\t1\t1\t+\t0\t100\t0.95\tfirst\n"
    )

    df_atoms = read_geese(geese_file=geese_file)

    assert list(df_atoms.columns) == ["name", "atom_nr", "class", "strand", "start", "end", "score", "note"]
    assert df_atoms.loc[0, "atom_nr"] == 1
    assert df_atoms.loc[0, "strand"] == "+"
    assert df_atoms.loc[0, "score"] == pytest.approx(0.95)
    assert df_atoms.loc[0, "note"] == "first"


# ---------------------------------------------------------------------
# Test: specific missing required columns are reported
# ---------------------------------------------------------------------
@pytest.mark.parametrize(
    ("file_content", "error_message"),
    [
        (
            "#name\tstart\tend\n"
            "gene1\t0\t100\n",
            "Missing required columns: class",
        ),
        (
            "#name\tclass\tend\n"
            "gene1\t1\t100\n",
            "Missing required columns: start",
        ),
        (
            "#name\tclass\tstart\n"
            "gene1\t1\t0\n",
            "Missing required columns: end",
        ),
        (
            "class\tstart\tend\n"
            "1\t0\t100\n",
            "Missing required columns: name",
        ),
        (
            "#name\tclass\n"
            "gene1\t1\n",
            "Missing required columns: start, end",
        ),
    ],
)
def test_read_geese_missing_columns(tmp_path: Path, file_content: str, error_message: str):
    """read_geese should raise ValueError with exact missing required column names."""
    broken_file = tmp_path / "broken.geese"
    broken_file.write_text(file_content)

    with pytest.raises(ValueError, match=error_message):
        read_geese(geese_file=broken_file)


# ---------------------------------------------------------------------
# Test: empty GEESE file is rejected
# ---------------------------------------------------------------------
def test_read_geese_empty_file(tmp_path: Path):
    """read_geese should raise ValueError if the GEESE file is empty."""
    empty_file = tmp_path / "empty.geese"
    empty_file.write_text("")

    with pytest.raises(ValueError, match=f"Malformed GEESE file: {empty_file}"):
        read_geese(geese_file=empty_file)


# ---------------------------------------------------------------------
# Test: malformed GEESE file is rejected
# ---------------------------------------------------------------------
def test_read_geese_malformed_file(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    """read_geese should raise ValueError if pandas cannot parse the GEESE file."""
    geese_file = tmp_path / "malformed.geese"
    geese_file.write_text("#name\tclass\tstart\tend\n")

    def mock_read_csv(*_args, **_kwargs):
        raise pd.errors.ParserError("broken parser state")

    monkeypatch.setattr("atomization_scorer.data_processing.geese_reader.pd.read_csv", mock_read_csv)

    with pytest.raises(ValueError, match=f"Malformed GEESE file: {geese_file}"):
        read_geese(geese_file=geese_file)


# ---------------------------------------------------------------------
# Test: missing GEESE file is rejected
# ---------------------------------------------------------------------
def test_read_geese_missing_file(tmp_path: Path):
    """read_geese should raise FileNotFoundError if the GEESE file does not exist."""
    missing_file = tmp_path / "missing.geese"

    with pytest.raises(FileNotFoundError):
        read_geese(geese_file=missing_file)
