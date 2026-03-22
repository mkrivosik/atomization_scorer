"""
geese_reader.py

Provides a function to load GEESE atomization data from a tabular file
into a pandas DataFrame, with validation of required columns.

Functions
---------
read_geese : Load the GEESE atomization tabular file into a pandas DataFrame
             with column renaming and required column check.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
from __future__ import annotations
from pathlib import Path

import pandas as pd
from pandas.errors import EmptyDataError, ParserError

from .utils import (
    check_required_columns,
    rename_column,
)

# --------------------------------------------------------------------------------------
# GEESE Reader
# --------------------------------------------------------------------------------------

def read_geese(geese_file: Path) -> pd.DataFrame:
    """
    Load a GEESE tabular file into a pandas DataFrame.

    The input may use either '#name' or 'name'; both are normalized to 'name'.
    Checks for required columns after normalization.

    Parameters
    ----------
    geese_file : Path
        Path to the input GEESE tabular file.

    Raises
    ------
    FileNotFoundError
        Raised if the GEESE file does not exist.
    ValueError
        Raised if the GEESE file is malformed or any required column is missing.

    Returns
    -------
    pd.DataFrame
        GEESE atomization data with required columns:
        ['name', 'class', 'start', 'end'].
        Additional columns from the input file are preserved.
    """
    if not geese_file.is_file():
        raise FileNotFoundError(f"GEESE file not found: {geese_file}")

    try:
        df = pd.read_csv(geese_file, sep=r"\s+", engine="python", comment=None)
    except (EmptyDataError, ParserError, UnicodeDecodeError) as error:
        raise ValueError(f"Malformed GEESE file: {geese_file}") from error

    df = rename_column(df=df, old_name="#name", new_name="name")

    required_columns = ["name", "class", "start", "end"]
    check_required_columns(df=df, required_columns=required_columns)

    return df
