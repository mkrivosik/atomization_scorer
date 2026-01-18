"""
geese_reader.py

Provides a function to load GEESE atomization data from a TSV file
into a pandas DataFrame, with validation of required columns.

Functions
---------
read_geese : Load GEESE atomization TSV file into a pandas DataFrame
             with column renaming and required column check.
"""

# ---------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------
import pandas as pd
from pathlib import Path
from .utils import (
    check_required_columns,
    rename_column,
)

# ---------------------------------------------------------------------
# GEESE Reader
# ---------------------------------------------------------------------

def read_geese(geese_file: Path) -> pd.DataFrame:
    """
    Loads a GEESE TSV file into a pandas DataFrame.

    Renames a '#name' column to 'name' for consistency and checks for required columns.

    Parameters
    ----------
    geese_file : Path
        Path to the GEESE TSV file.

    Raises
    ------
    FileNotFoundError
        Raised if the GEESE file does not exist.
    ValueError
        Raised if any required column is missing in the input file.

    Returns
    -------
    pd.DataFrame
        GEESE atomization data with columns:
        ['name', 'atom_nr', 'class', 'strand', 'start', 'end']
    """
    if not geese_file.is_file():
        raise FileNotFoundError(f"GEESE file not found: {geese_file}")

    df = pd.read_csv(geese_file, sep=r"\s+", engine="python", comment=None)
    df = rename_column(df=df, old_name="#name", new_name="name")

    required_columns = ["name", "class", "start", "end"]
    check_required_columns(df=df, required_columns=required_columns)

    return df
