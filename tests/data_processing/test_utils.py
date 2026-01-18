"""
Tests for the utils.py functions: check_required_columns and rename_column.
"""

import pytest
import pandas as pd
from atomization_scorer.data_processing import check_required_columns, rename_column


# ---------------------------------------------------------------------
# Test: check_required_columns passes if all required columns exist
# ---------------------------------------------------------------------
def test_check_required_columns_valid():
    """check_required_columns should not raise an error if all required columns exist."""
    df = pd.DataFrame({
        'A': [1, 2],
        'B': [3, 4],
        'C': [5, 6]
    })
    required = ['A', 'B']

    check_required_columns(df=df, required_columns=required)


# ---------------------------------------------------------------------
# Test: check_required_columns raises ValueError if columns are missing
# ---------------------------------------------------------------------
def test_check_required_columns_missing():
    """check_required_columns should raise ValueError if required columns are missing."""
    df = pd.DataFrame({
        'A': [1, 2],
        'B': [3, 4]
    })
    required = ['A', 'B', 'C']

    with pytest.raises(ValueError) as exception_info:
        check_required_columns(df=df, required_columns=required)

    assert "Missing required columns: C" in str(exception_info.value)


# ---------------------------------------------------------------------
# Test: rename_column renames a column if it exists
# ---------------------------------------------------------------------
def test_rename_column_existing():
    """rename_column should rename the column if it exists in the DataFrame."""
    df = pd.DataFrame({
        'old_name': [1, 2],
        'other': [3, 4]
    })
    df = rename_column(df=df, old_name='old_name', new_name='new_name')

    assert 'new_name' in df.columns
    assert 'old_name' not in df.columns
    assert df['new_name'].tolist() == [1, 2]


# ---------------------------------------------------------------------
# Test: rename_column with empty DataFrame
# ---------------------------------------------------------------------
def test_rename_column_empty_dataframe():
    """rename_column should not fail on an empty DataFrame and return it unchanged."""
    df = pd.DataFrame()
    df_renamed = rename_column(df=df, old_name='any_column', new_name='new_name')

    assert df_renamed.empty
    assert df_renamed.columns.tolist() == []


# ---------------------------------------------------------------------
# Test: rename_column does nothing if column does not exist
# ---------------------------------------------------------------------
def test_rename_column_nonexistent():
    """rename_column should not modify DataFrame if old_name does not exist."""
    df = pd.DataFrame({
        'A': [1, 2],
        'B': [3, 4]
    })
    df_copy = df.copy()
    df = rename_column(df=df, old_name='C', new_name='D')  # 'C' does not exist

    assert df.equals(df_copy)
