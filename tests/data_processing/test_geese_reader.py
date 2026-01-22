"""
Tests for the read_geese() function.
"""

from pathlib import Path

import pytest

from atomization_scorer import read_geese


# ---------------------------------------------------------------------
# Test: loads a valid GEESE TSV and contains all required columns
# ---------------------------------------------------------------------
def test_read_geese_valid(tmp_path: Path):
    """read_geese should load a valid GEESE TSV and contain all required columns."""
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
# Test: raises ValueError if required columns are missing
# ---------------------------------------------------------------------
def test_read_geese_missing_columns(tmp_path: Path):
    """read_geese should raise ValueError if required columns are missing."""
    broken_file = tmp_path / "broken.geese"
    broken_file.write_text("apple\torange\tpear\n")  # invalid GEESE format

    with pytest.raises(ValueError):
        read_geese(geese_file=broken_file)


# ---------------------------------------------------------------------
# Test: handles empty GEESE file
# ---------------------------------------------------------------------
def test_read_geese_empty_file(tmp_path: Path):
    """read_geese should raise ValueError if the GEESE file is empty."""
    empty_file = tmp_path / "empty.geese"
    empty_file.write_text("")

    with pytest.raises(ValueError):
        read_geese(geese_file=empty_file)


# ---------------------------------------------------------------------
# Test: raises FileNotFoundError if file does not exist
# ---------------------------------------------------------------------
def test_read_geese_missing_file(tmp_path: Path):
    """read_geese should raise FileNotFoundError if the file does not exist."""
    missing_file = tmp_path / "missing.geese"

    with pytest.raises(FileNotFoundError):
        read_geese(geese_file=missing_file)
