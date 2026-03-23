"""
utils.py

Shared utility functions for data processing in genomes atomization scoring.

Functions
---------
check_required_columns  : Check that all required columns are present in a DataFrame.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
from __future__ import annotations
import pandas as pd

# --------------------------------------------------------------------------------------
# Utils Functions
# --------------------------------------------------------------------------------------

def check_required_columns(df: pd.DataFrame, required_columns: list[str]) -> None:
    """
    Check that all required columns are present in a DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Table to check.
    required_columns : list of str
        Columns that must exist.

    Raises
    ------
    ValueError
        Raised if any required column is missing.

    Returns
    -------
    None
    """
    missing = [column for column in required_columns if column not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")
