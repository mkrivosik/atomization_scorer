"""
Tests for the artificial_degradation.py module.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
from pathlib import Path

import pandas as pd
import pytest

from atomization_scorer import (
    compute_overall_score,
    compute_true_alignment,
    degrade_atomization,
    read_geese,
    write_geese,
)


# --------------------------------------------------------------------------------------
# Helper: build a valid atomization table for degradation tests
# --------------------------------------------------------------------------------------
def _build_atomization_dataframe() -> pd.DataFrame:
    """Create a valid GEESE DataFrame with several distinct classes."""
    return pd.DataFrame(
        {
            "name": ["genome1", "genome1", "genome1", "genome1", "genome1"],
            "class": [1, 2, 3, 1, 2],
            "start": [0, 10, 20, 30, 40],
            "end": [10, 20, 30, 40, 50],
        }
    )


# --------------------------------------------------------------------------------------
# Test: degrade_atomization rejects invalid degradation fractions
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    "degradation_fraction",
    [-0.1, 1.1],
)
def test_degrade_atomization_invalid_fraction(tmp_path: Path, degradation_fraction: float):
    """degrade_atomization should reject degradation fractions outside [0.0, 1.0]."""
    atomization_file = tmp_path / "input.geese"
    write_geese(df_atoms=_build_atomization_dataframe(), output_path=atomization_file)

    with pytest.raises(ValueError, match="degradation_fraction must be between 0.0 and 1.0"):
        degrade_atomization(
            atomization_file=atomization_file,
            output_directory=tmp_path,
            degradation_fraction=degradation_fraction,
            random_seed=7,
        )


# --------------------------------------------------------------------------------------
# Test: degrade_atomization rejects atomization with only one distinct class
# --------------------------------------------------------------------------------------
def test_degrade_atomization_single_class_rejected(tmp_path: Path):
    """degrade_atomization should raise ValueError if there is no alternative class to assign."""
    atomization_file = tmp_path / "input.geese"
    df_atoms = pd.DataFrame(
        {
            "name": ["genome1", "genome1"],
            "class": [1, 1],
            "start": [0, 10],
            "end": [10, 20],
        }
    )
    write_geese(df_atoms=df_atoms, output_path=atomization_file)

    with pytest.raises(ValueError, match="requires at least two distinct atomization classes"):
        degrade_atomization(
            atomization_file=atomization_file,
            output_directory=tmp_path,
            degradation_fraction=0.5,
            random_seed=7,
        )


# --------------------------------------------------------------------------------------
# Test: zero degradation preserves the original atomization exactly
# --------------------------------------------------------------------------------------
def test_degrade_atomization_zero_fraction_round_trip(tmp_path: Path):
    """degrade_atomization should preserve all rows when degradation_fraction is zero."""
    atomization_file = tmp_path / "input.geese"
    df_atoms = _build_atomization_dataframe()
    write_geese(df_atoms=df_atoms, output_path=atomization_file)

    path = degrade_atomization(
        atomization_file=atomization_file,
        output_directory=tmp_path,
        degradation_fraction=0.0,
        random_seed=7,
    )

    result = read_geese(geese_file=path)

    assert path == tmp_path / "degraded_atomization.geese"
    pd.testing.assert_frame_equal(result, df_atoms)


# --------------------------------------------------------------------------------------
# Test: degrade_atomization changes classes and preserves coordinates
# --------------------------------------------------------------------------------------
def test_degrade_atomization_changes_selected_classes_only(tmp_path: Path):
    """degrade_atomization should change only class values for the selected atom fraction."""
    atomization_file = tmp_path / "input.geese"
    df_atoms = _build_atomization_dataframe()
    write_geese(df_atoms=df_atoms, output_path=atomization_file)

    output_file = degrade_atomization(
        atomization_file=atomization_file,
        output_directory=tmp_path,
        degradation_fraction=0.4,
        random_seed=7,
    )

    result = read_geese(geese_file=output_file)
    changed_rows = result["class"] != df_atoms["class"]

    assert changed_rows.sum() == 2
    assert result.loc[changed_rows, "class"].tolist() != df_atoms.loc[changed_rows, "class"].tolist()
    pd.testing.assert_series_equal(result["name"], df_atoms["name"], check_names=True)
    pd.testing.assert_series_equal(result["start"], df_atoms["start"], check_names=True)
    pd.testing.assert_series_equal(result["end"], df_atoms["end"], check_names=True)


# --------------------------------------------------------------------------------------
# Test: degrade_atomization is reproducible with a fixed random seed
# --------------------------------------------------------------------------------------
def test_degrade_atomization_is_reproducible_with_seed(tmp_path: Path):
    """degrade_atomization should produce the same degraded classes for the same random seed."""
    atomization_file = tmp_path / "input.geese"
    first_output_directory = tmp_path / "first_output"
    second_output_directory = tmp_path / "second_output"
    write_geese(df_atoms=_build_atomization_dataframe(), output_path=atomization_file)

    first_output = degrade_atomization(
        atomization_file=atomization_file,
        output_directory=first_output_directory,
        degradation_fraction=0.6,
        random_seed=11,
    )
    second_output = degrade_atomization(
        atomization_file=atomization_file,
        output_directory=second_output_directory,
        degradation_fraction=0.6,
        random_seed=11,
    )

    first_result = read_geese(geese_file=first_output)
    second_result = read_geese(geese_file=second_output)

    pd.testing.assert_frame_equal(first_result, second_result)


# --------------------------------------------------------------------------------------
# Test: overall score does not increase as degradation grows
# --------------------------------------------------------------------------------------
def test_degrade_atomization_real_pipeline_overall_score_non_increasing(
    mini_fasta: Path,
    mini_geese: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    """compute_overall_score should not increase when class degradation is increased."""
    degradation_fractions = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    reference_true_geese = compute_true_alignment(
        genomes_file=mini_fasta,
        atomization_file=mini_geese,
        output_directory=tmp_path / "reference_true_alignment",
    )
    monkeypatch.setattr(
        "atomization_scorer.scoring_system.alignment_score.compute_true_alignment",
        lambda **_kwargs: reference_true_geese,
    )
    previous_score = None

    for degradation_fraction in degradation_fractions:
        degraded_file = degrade_atomization(
            atomization_file=mini_geese,
            output_directory=tmp_path / f"degraded_{degradation_fraction:.1f}",
            degradation_fraction=degradation_fraction,
            random_seed=7,
        )
        current_score = compute_overall_score(
            genomes_file=mini_fasta,
            atomization_file=degraded_file,
            output_directory=tmp_path / f"score_{degradation_fraction:.1f}",
        )

        if previous_score is not None:
            assert current_score <= previous_score

        previous_score = current_score
