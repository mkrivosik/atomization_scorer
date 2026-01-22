"""
Tests for the compute_true_alignment() function.
"""

from pathlib import Path
from unittest.mock import patch

import pytest

from atomization_scorer import compute_true_alignment


# -----------------------------------------------------------------------------
# Test: calls all pipeline steps and returns GEESE file
# -----------------------------------------------------------------------------
@patch("atomization_scorer.pipeline.true_pipeline.extract_representatives")
@patch("atomization_scorer.pipeline.true_pipeline.align_with_minimap2")
@patch("atomization_scorer.pipeline.true_pipeline.filter_paf")
@patch("atomization_scorer.pipeline.true_pipeline.paf_to_geese")
def test_compute_true_alignment_pipeline(
    mock_paf_to_geese,
    mock_filter_paf,
    mock_minimap2,
    mock_extract,
    mini_fasta: Path,
    mini_geese: Path,
    output_dir: Path
):
    """Test that compute_true_alignment calls all pipeline steps with correct arguments."""
    # Run pipeline
    geese_path = compute_true_alignment(
        genomes_file=mini_fasta,
        atomization_file=mini_geese,
        output_directory=output_dir,
        mode="mash",
        minimum_similarity=0.95,
        minimum_alignment_length=500
    )

    # Expected paths
    representatives_fasta = output_dir / "mash_representatives.fa"
    paf_file = output_dir / "minimap2_alignments.paf"
    filtered_paf = output_dir / "minimap2_alignment_filtered.paf"
    true_geese = output_dir / "true_atomization.geese"

    assert geese_path == true_geese

    # Check calls
    mock_extract.assert_called_once_with(
        genomes_file=mini_fasta,
        atomization_file=mini_geese,
        output_path=representatives_fasta,
        mode="mash"
    )

    mock_minimap2.assert_called_once_with(
        query=mini_fasta,
        target=representatives_fasta,
        output_path=paf_file
    )

    mock_filter_paf.assert_called_once_with(
        paf_file=paf_file,
        output_file=filtered_paf,
        minimum_similarity=0.95,
        minimum_alignment_length=500
    )

    mock_paf_to_geese.assert_called_once_with(
        paf_file=filtered_paf,
        output_file=true_geese
    )


# -----------------------------------------------------------------------------
# Test: works with "first" mode instead of "mash"
# -----------------------------------------------------------------------------
@patch("atomization_scorer.pipeline.true_pipeline.extract_representatives")
@patch("atomization_scorer.pipeline.true_pipeline.align_with_minimap2")
@patch("atomization_scorer.pipeline.true_pipeline.filter_paf")
@patch("atomization_scorer.pipeline.true_pipeline.paf_to_geese")
def test_compute_true_alignment_first_mode(
    mock_paf_to_geese,
    mock_filter_paf,
    mock_minimap2,
    mock_extract,
    mini_fasta: Path,
    mini_geese: Path,
    output_dir: Path
):
    """Test that compute_true_alignment works when mode is set to 'first'."""
    # Run pipeline
    geese_path = compute_true_alignment(
        genomes_file=mini_fasta,
        atomization_file=mini_geese,
        output_directory=output_dir,
        mode="first",
        minimum_similarity=0.95,
        minimum_alignment_length=500
    )

    representatives_fasta = output_dir / "first_representatives.fa"
    true_geese = output_dir / "true_atomization.geese"

    assert geese_path == true_geese

    mock_extract.assert_called_once_with(
        genomes_file=mini_fasta,
        atomization_file=mini_geese,
        output_path=representatives_fasta,
        mode="first"
    )

    mock_minimap2.assert_called_once()
    mock_filter_paf.assert_called_once()
    mock_paf_to_geese.assert_called_once()


# -----------------------------------------------------------------------------
# Test: raises FileNotFoundError if genomes file is missing
# -----------------------------------------------------------------------------
def test_compute_true_alignment_missing_genomes(mini_geese: Path, tmp_path: Path, output_dir: Path):
    """compute_true_alignment should raise FileNotFoundError if the genomes FASTA file does not exist."""
    missing_genomes = tmp_path / "missing_genomes.fa"

    with pytest.raises(FileNotFoundError, match="Genomes FASTA file not found"):
        compute_true_alignment(
            genomes_file=missing_genomes,
            atomization_file=mini_geese,
            output_directory=output_dir
        )


# -----------------------------------------------------------------------------
# Test: raises FileNotFoundError if atomization file is missing
# -----------------------------------------------------------------------------
def test_compute_true_alignment_missing_atomization(mini_fasta: Path, tmp_path: Path, output_dir: Path):
    """compute_true_alignment should raise FileNotFoundError if the atomization file does not exist."""
    missing_atomization = tmp_path / "missing_atomization.geese"

    with pytest.raises(FileNotFoundError, match="Atomization file not found"):
        compute_true_alignment(
            genomes_file=mini_fasta,
            atomization_file=missing_atomization,
            output_directory=output_dir
        )
