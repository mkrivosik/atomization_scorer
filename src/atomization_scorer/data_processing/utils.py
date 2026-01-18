"""
utils.py

Shared utility functions for data processing in genomes atomization scoring.

Functions
---------
check_required_columns  : Check that all required columns are present in a DataFrame.
rename_column           : Rename a column in a DataFrame if it exists.
"""

# ---------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------
import pandas as pd

# ---------------------------------------------------------------------
# Utils functions
# ---------------------------------------------------------------------

def check_required_columns(df: pd.DataFrame, required_columns: list[str]) -> None:
    """
    Check that all required columns are present in a DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to check.
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


def rename_column(df: pd.DataFrame, old_name: str, new_name: str) -> pd.DataFrame:
    """
    Rename a column if it exists.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to modify.
    old_name : str
        Column name to rename.
    new_name : str
        New column name.

    Returns
    -------
    pd.DataFrame
        Modified DataFrame.
    """
    if old_name in df.columns:
        df.rename(columns={old_name: new_name}, inplace=True)
    return df
