"""
Tests for the compute_true_alignment() function.
"""

from pathlib import Path
from unittest.mock import patch

import pytest

from atomization_scorer import compute_true_alignment


# --------------------------------------------------------------------------------------
# Test: calls all pipeline steps and returns GEESE file
# --------------------------------------------------------------------------------------
@patch("atomization_scorer.pipeline.true_pipeline.plot_genome_atomization")
@patch("atomization_scorer.pipeline.true_pipeline.paf_to_geese")
@patch("atomization_scorer.pipeline.true_pipeline.filter_paf")
@patch("atomization_scorer.pipeline.true_pipeline.align_with_minimap2")
@patch("atomization_scorer.pipeline.true_pipeline.extract_representatives")
def test_compute_true_alignment_pipeline(
    mock_extract,
    mock_minimap2,
    mock_filter_paf,
    mock_paf_to_geese,
    mock_visualization,
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
    visualization_directory = output_dir / "atomization_visualization"

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

    mock_visualization.assert_called_once_with(
        genomes_file=mini_fasta,
        true_atoms_file=true_geese,
        predicted_atoms_file=mini_geese,
        output_directory=visualization_directory
    )


# --------------------------------------------------------------------------------------
# Test: works with "first" mode instead of "mash"
# --------------------------------------------------------------------------------------
@patch("atomization_scorer.pipeline.true_pipeline.plot_genome_atomization")
@patch("atomization_scorer.pipeline.true_pipeline.paf_to_geese")
@patch("atomization_scorer.pipeline.true_pipeline.filter_paf")
@patch("atomization_scorer.pipeline.true_pipeline.align_with_minimap2")
@patch("atomization_scorer.pipeline.true_pipeline.extract_representatives")
def test_compute_true_alignment_first_mode(
    mock_extract,
    mock_minimap2,
    mock_filter_paf,
    mock_paf_to_geese,
    mock_visualization,
    mini_fasta: Path,
    mini_geese: Path,
    output_dir: Path
):
    """Test that compute_true_alignment works when the mode is set to 'first'."""
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
    mock_visualization.assert_called_once()


# --------------------------------------------------------------------------------------
# Test: forwards non-default filter parameters
# --------------------------------------------------------------------------------------
@patch("atomization_scorer.pipeline.true_pipeline.plot_genome_atomization")
@patch("atomization_scorer.pipeline.true_pipeline.paf_to_geese")
@patch("atomization_scorer.pipeline.true_pipeline.filter_paf")
@patch("atomization_scorer.pipeline.true_pipeline.align_with_minimap2")
@patch("atomization_scorer.pipeline.true_pipeline.extract_representatives")
def test_compute_true_alignment_forwards_custom_filter_parameters(
    mock_extract,
    mock_minimap2,
    mock_filter_paf,
    mock_paf_to_geese,
    mock_visualization,
    mini_fasta: Path,
    mini_geese: Path,
    output_dir: Path
):
    """compute_true_alignment should forward non-default filter parameters to PAF filtering."""
    compute_true_alignment(
        genomes_file=mini_fasta,
        atomization_file=mini_geese,
        output_directory=output_dir,
        minimum_similarity=0.8,
        minimum_alignment_length=1234
    )

    filtered_paf = output_dir / "minimap2_alignment_filtered.paf"
    paf_file = output_dir / "minimap2_alignments.paf"

    mock_extract.assert_called_once()
    mock_minimap2.assert_called_once()
    mock_filter_paf.assert_called_once_with(
        paf_file=paf_file,
        output_file=filtered_paf,
        minimum_similarity=0.8,
        minimum_alignment_length=1234
    )
    mock_paf_to_geese.assert_called_once()
    mock_visualization.assert_called_once()


# --------------------------------------------------------------------------------------
# Test: creates output directory if missing
# --------------------------------------------------------------------------------------
@patch("atomization_scorer.pipeline.true_pipeline.plot_genome_atomization")
@patch("atomization_scorer.pipeline.true_pipeline.paf_to_geese")
@patch("atomization_scorer.pipeline.true_pipeline.filter_paf")
@patch("atomization_scorer.pipeline.true_pipeline.align_with_minimap2")
@patch("atomization_scorer.pipeline.true_pipeline.extract_representatives")
def test_compute_true_alignment_creates_output_directory(
    mock_extract,
    mock_minimap2,
    mock_filter_paf,
    mock_paf_to_geese,
    mock_visualization,
    mini_fasta: Path,
    mini_geese: Path,
    tmp_path: Path
):
    """compute_true_alignment should create a missing output directory before running the pipeline."""
    output_dir = tmp_path / "nested" / "true_pipeline_output"

    assert not output_dir.exists()

    compute_true_alignment(
        genomes_file=mini_fasta,
        atomization_file=mini_geese,
        output_directory=output_dir
    )

    assert output_dir.is_dir()
    mock_extract.assert_called_once()
    mock_minimap2.assert_called_once()
    mock_filter_paf.assert_called_once()
    mock_paf_to_geese.assert_called_once()
    mock_visualization.assert_called_once()


# --------------------------------------------------------------------------------------
# Test: raises FileNotFoundError if an input file is missing
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("missing_side", "expected_message"),
    [
        ("genomes", "Genomes FASTA file not found"),
        ("atomization", "Atomization file not found"),
    ],
)
def test_compute_true_alignment_missing_input_file(
    mini_fasta: Path,
    mini_geese: Path,
    tmp_path: Path,
    output_dir: Path,
    missing_side: str,
    expected_message: str,
):
    """compute_true_alignment should raise FileNotFoundError if an input file does not exist."""
    genomes_file = mini_fasta
    atomization_file = mini_geese

    if missing_side == "genomes":
        genomes_file = tmp_path / "missing_genomes.fa"
    else:
        atomization_file = tmp_path / "missing_atomization.geese"

    with pytest.raises(FileNotFoundError, match=expected_message):
        compute_true_alignment(
            genomes_file=genomes_file,
            atomization_file=atomization_file,
            output_directory=output_dir
        )


# --------------------------------------------------------------------------------------
# Test: propagates minimap2 alignment failure
# --------------------------------------------------------------------------------------
@patch("atomization_scorer.pipeline.true_pipeline.plot_genome_atomization")
@patch("atomization_scorer.pipeline.true_pipeline.paf_to_geese")
@patch("atomization_scorer.pipeline.true_pipeline.filter_paf")
@patch("atomization_scorer.pipeline.true_pipeline.align_with_minimap2")
@patch("atomization_scorer.pipeline.true_pipeline.extract_representatives")
def test_compute_true_alignment_propagates_alignment_failure(
    mock_extract,
    mock_minimap2,
    mock_filter_paf,
    mock_paf_to_geese,
    mock_visualization,
    mini_fasta: Path,
    mini_geese: Path,
    output_dir: Path
):
    """compute_true_alignment should propagate minimap2 failures and stop later stages."""
    mock_minimap2.side_effect = RuntimeError("alignment failed")

    with pytest.raises(RuntimeError, match="alignment failed"):
        compute_true_alignment(
            genomes_file=mini_fasta,
            atomization_file=mini_geese,
            output_directory=output_dir
        )

    mock_extract.assert_called_once()
    mock_minimap2.assert_called_once()
    mock_filter_paf.assert_not_called()
    mock_paf_to_geese.assert_not_called()
    mock_visualization.assert_not_called()


# --------------------------------------------------------------------------------------
# Test: propagates PAF filtering failure
# --------------------------------------------------------------------------------------
@patch("atomization_scorer.pipeline.true_pipeline.plot_genome_atomization")
@patch("atomization_scorer.pipeline.true_pipeline.paf_to_geese")
@patch("atomization_scorer.pipeline.true_pipeline.filter_paf")
@patch("atomization_scorer.pipeline.true_pipeline.align_with_minimap2")
@patch("atomization_scorer.pipeline.true_pipeline.extract_representatives")
def test_compute_true_alignment_propagates_filter_failure(
    mock_extract,
    mock_minimap2,
    mock_filter_paf,
    mock_paf_to_geese,
    mock_visualization,
    mini_fasta: Path,
    mini_geese: Path,
    output_dir: Path
):
    """compute_true_alignment should propagate PAF filtering failures and stop later stages."""
    mock_filter_paf.side_effect = RuntimeError("filter failed")

    with pytest.raises(RuntimeError, match="filter failed"):
        compute_true_alignment(
            genomes_file=mini_fasta,
            atomization_file=mini_geese,
            output_directory=output_dir
        )

    mock_extract.assert_called_once()
    mock_minimap2.assert_called_once()
    mock_filter_paf.assert_called_once()
    mock_paf_to_geese.assert_not_called()
    mock_visualization.assert_not_called()


# --------------------------------------------------------------------------------------
# Test: propagates visualization failure
# --------------------------------------------------------------------------------------
@patch("atomization_scorer.pipeline.true_pipeline.plot_genome_atomization")
@patch("atomization_scorer.pipeline.true_pipeline.paf_to_geese")
@patch("atomization_scorer.pipeline.true_pipeline.filter_paf")
@patch("atomization_scorer.pipeline.true_pipeline.align_with_minimap2")
@patch("atomization_scorer.pipeline.true_pipeline.extract_representatives")
def test_compute_true_alignment_propagates_visualization_failure(
    mock_extract,
    mock_minimap2,
    mock_filter_paf,
    mock_paf_to_geese,
    mock_visualization,
    mini_fasta: Path,
    mini_geese: Path,
    output_dir: Path
):
    """compute_true_alignment should propagate visualization failures after earlier steps succeed."""
    mock_visualization.side_effect = RuntimeError("visualization failed")

    with pytest.raises(RuntimeError, match="visualization failed"):
        compute_true_alignment(
            genomes_file=mini_fasta,
            atomization_file=mini_geese,
            output_directory=output_dir
        )

    mock_extract.assert_called_once()
    mock_minimap2.assert_called_once()
    mock_filter_paf.assert_called_once()
    mock_paf_to_geese.assert_called_once()
    mock_visualization.assert_called_once()
