"""
geese_writer.py

Provides a function to write GEESE atomization data from a pandas DataFrame
into a tabular file, with validation of required columns.

Functions
---------
write_geese : Write GEESE atomization tabular data to a file while preserving
              additional columns.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
from __future__ import annotations
from pathlib import Path

import pandas as pd

from .utils import check_required_columns


# --------------------------------------------------------------------------------------
# GEESE Writer
# --------------------------------------------------------------------------------------
def write_geese(df_atoms: pd.DataFrame, output_path: Path) -> Path:
    """
    Write GEESE atomization data from a pandas DataFrame into a tabular file.

    Parameters
    ----------
    df_atoms : pd.DataFrame
        GEESE atomization data containing required columns ['name', 'class',
        'start', 'end']. Additional columns are preserved.
    output_path : Path
        Path where the GEESE file should be written.

    Raises
    ------
    ValueError
        Raised if any required GEESE column is missing.

    Returns
    -------
    Path
        Written GEESE file path.
    """
    required_columns = ["name", "class", "start", "end"]
    check_required_columns(df=df_atoms, required_columns=required_columns)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    output_df = df_atoms.copy()
    output_df = output_df.rename(columns={"name": "#name"})
    output_df.to_csv(output_path, sep="\t", index=False)

    return output_path
