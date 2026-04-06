"""
utils.py

Shared utility functions for data processing in genomes atomization scoring.

Functions
---------
check_required_columns  : Check that all required columns are present in a DataFrame.
sanitize_path_component : Convert an arbitrary identifier into a filesystem-safe name.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
from __future__ import annotations
import hashlib
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


def sanitize_path_component(value: str, fallback_prefix: str = "component") -> str:
    """
    Convert an arbitrary identifier into a filesystem-safe name.

    Parameters
    ----------
    value : str
        Input identifier to sanitize.
    fallback_prefix : str, optional, default="component"
        Prefix used when the sanitized result is empty.

    Returns
    -------
    str
        Filesystem-safe name.
    """
    safe_characters = {"-", "_", "."}
    sanitized = "".join(
        character if character.isalnum() or character in safe_characters else "_"
        for character in value
    ).strip("_")

    if sanitized:
        return sanitized

    digest = hashlib.sha1(value.encode("utf-8")).hexdigest()[:8]
    return f"{fallback_prefix}_{digest}"
