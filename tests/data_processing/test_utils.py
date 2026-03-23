"""
Tests for the utils.py functions.
"""

import pandas as pd
import pytest

from atomization_scorer.data_processing import check_required_columns


# --------------------------------------------------------------------------------------
# Test: check_required_columns passes if all required columns exist
# --------------------------------------------------------------------------------------
def test_check_required_columns_valid():
    """check_required_columns should not raise an error if all required columns exist."""
    df = pd.DataFrame({
        'A': [1, 2],
        'B': [3, 4],
        'C': [5, 6]
    })
    required = ['A', 'B']

    check_required_columns(df=df, required_columns=required)


# --------------------------------------------------------------------------------------
# Test: check_required_columns raises ValueError if one column is missing
# --------------------------------------------------------------------------------------
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


# --------------------------------------------------------------------------------------
# Test: check_required_columns raises ValueError if several columns are missing
# --------------------------------------------------------------------------------------
def test_check_required_columns_multiple_missing():
    """check_required_columns should raise ValueError if several required columns are missing."""
    df = pd.DataFrame({
        'A': [1, 2],
    })
    required = ['A', 'B', 'C', 'D']

    with pytest.raises(ValueError) as exception_info:
        check_required_columns(df=df, required_columns=required)

    assert "Missing required columns: B, C, D" in str(exception_info.value)
