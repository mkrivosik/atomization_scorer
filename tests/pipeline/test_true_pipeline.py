"""
Tests for the true_pipeline.py module.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
from pathlib import Path
from unittest.mock import patch

import pytest

from atomization_scorer import compute_true_alignment


# --------------------------------------------------------------------------------------
# Test: calls all pipeline steps and returns GEESE file
# --------------------------------------------------------------------------------------
@patch("atomization_scorer.pipeline.true_pipeline.paf_to_geese")
@patch("atomization_scorer.pipeline.true_pipeline.validate_non_overlapping_geese")
@patch("atomization_scorer.pipeline.true_pipeline.resolve_paf_overlaps")
@patch("atomization_scorer.pipeline.true_pipeline.diagnose_paf_overlaps")
@patch("atomization_scorer.pipeline.true_pipeline.filter_paf")
@patch("atomization_scorer.pipeline.true_pipeline.align_with_minimap2")
@patch("atomization_scorer.pipeline.true_pipeline.extract_representatives")
def test_compute_true_alignment_pipeline(
    mock_extract,
    mock_minimap2,
    mock_filter_paf,
    mock_diagnose,
    mock_resolve_paf,
    mock_validate_geese,
    mock_paf_to_geese,
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
        minimum_alignment_length=500,
        run_overlap_diagnostics=False,
    )

    # Expected paths
    representatives_fasta = output_dir / "mash_representatives.fa"
    paf_file = output_dir / "minimap2_alignments.paf"
    filtered_paf = output_dir / "minimap2_alignment_filtered.paf"
    resolved_paf = output_dir / "minimap2_alignment_resolved.paf"
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
    mock_diagnose.assert_not_called()

    mock_resolve_paf.assert_called_once_with(
        paf_file=filtered_paf,
        output_file=resolved_paf,
        minimum_alignment_length=500,
    )

    mock_paf_to_geese.assert_called_once_with(
        paf_file=resolved_paf,
        output_file=true_geese
    )

    mock_validate_geese.assert_called_once_with(true_geese)


# --------------------------------------------------------------------------------------
# Test: works with "first" mode instead of "mash"
# --------------------------------------------------------------------------------------
@patch("atomization_scorer.pipeline.true_pipeline.paf_to_geese")
@patch("atomization_scorer.pipeline.true_pipeline.validate_non_overlapping_geese")
@patch("atomization_scorer.pipeline.true_pipeline.resolve_paf_overlaps")
@patch("atomization_scorer.pipeline.true_pipeline.diagnose_paf_overlaps")
@patch("atomization_scorer.pipeline.true_pipeline.filter_paf")
@patch("atomization_scorer.pipeline.true_pipeline.align_with_minimap2")
@patch("atomization_scorer.pipeline.true_pipeline.extract_representatives")
def test_compute_true_alignment_first_mode(
    mock_extract,
    mock_minimap2,
    mock_filter_paf,
    mock_diagnose,
    mock_resolve_paf,
    mock_validate_geese,
    mock_paf_to_geese,
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
        minimum_alignment_length=500,
        run_overlap_diagnostics=False,
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
    mock_diagnose.assert_not_called()
    mock_resolve_paf.assert_called_once()
    mock_paf_to_geese.assert_called_once()
    mock_validate_geese.assert_called_once()


# --------------------------------------------------------------------------------------
# Test: forwards non-default filter parameters
# --------------------------------------------------------------------------------------
@patch("atomization_scorer.pipeline.true_pipeline.paf_to_geese")
@patch("atomization_scorer.pipeline.true_pipeline.validate_non_overlapping_geese")
@patch("atomization_scorer.pipeline.true_pipeline.resolve_paf_overlaps")
@patch("atomization_scorer.pipeline.true_pipeline.diagnose_paf_overlaps")
@patch("atomization_scorer.pipeline.true_pipeline.filter_paf")
@patch("atomization_scorer.pipeline.true_pipeline.align_with_minimap2")
@patch("atomization_scorer.pipeline.true_pipeline.extract_representatives")
def test_compute_true_alignment_forwards_custom_filter_parameters(
    mock_extract,
    mock_minimap2,
    mock_filter_paf,
    mock_diagnose,
    mock_resolve_paf,
    mock_validate_geese,
    mock_paf_to_geese,
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
        minimum_alignment_length=1234,
        run_overlap_diagnostics=False,
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
    mock_diagnose.assert_not_called()
    mock_resolve_paf.assert_called_once()
    mock_paf_to_geese.assert_called_once()
    mock_validate_geese.assert_called_once()


# --------------------------------------------------------------------------------------
# Test: creates output directory if missing
# --------------------------------------------------------------------------------------
@patch("atomization_scorer.pipeline.true_pipeline.paf_to_geese")
@patch("atomization_scorer.pipeline.true_pipeline.validate_non_overlapping_geese")
@patch("atomization_scorer.pipeline.true_pipeline.resolve_paf_overlaps")
@patch("atomization_scorer.pipeline.true_pipeline.diagnose_paf_overlaps")
@patch("atomization_scorer.pipeline.true_pipeline.filter_paf")
@patch("atomization_scorer.pipeline.true_pipeline.align_with_minimap2")
@patch("atomization_scorer.pipeline.true_pipeline.extract_representatives")
def test_compute_true_alignment_creates_output_directory(
    mock_extract,
    mock_minimap2,
    mock_filter_paf,
    mock_diagnose,
    mock_resolve_paf,
    mock_validate_geese,
    mock_paf_to_geese,
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
        output_directory=output_dir,
        run_overlap_diagnostics=False,
    )

    assert output_dir.is_dir()
    mock_extract.assert_called_once()
    mock_minimap2.assert_called_once()
    mock_filter_paf.assert_called_once()
    mock_diagnose.assert_not_called()
    mock_resolve_paf.assert_called_once()
    mock_paf_to_geese.assert_called_once()
    mock_validate_geese.assert_called_once()


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
@patch("atomization_scorer.pipeline.true_pipeline.paf_to_geese")
@patch("atomization_scorer.pipeline.true_pipeline.validate_non_overlapping_geese")
@patch("atomization_scorer.pipeline.true_pipeline.resolve_paf_overlaps")
@patch("atomization_scorer.pipeline.true_pipeline.diagnose_paf_overlaps")
@patch("atomization_scorer.pipeline.true_pipeline.filter_paf")
@patch("atomization_scorer.pipeline.true_pipeline.align_with_minimap2")
@patch("atomization_scorer.pipeline.true_pipeline.extract_representatives")
def test_compute_true_alignment_propagates_alignment_failure(
    mock_extract,
    mock_minimap2,
    mock_filter_paf,
    mock_diagnose,
    mock_resolve_paf,
    mock_validate_geese,
    mock_paf_to_geese,
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
    mock_diagnose.assert_not_called()
    mock_resolve_paf.assert_not_called()
    mock_paf_to_geese.assert_not_called()
    mock_validate_geese.assert_not_called()


# --------------------------------------------------------------------------------------
# Test: propagates PAF filtering failure
# --------------------------------------------------------------------------------------
@patch("atomization_scorer.pipeline.true_pipeline.paf_to_geese")
@patch("atomization_scorer.pipeline.true_pipeline.validate_non_overlapping_geese")
@patch("atomization_scorer.pipeline.true_pipeline.resolve_paf_overlaps")
@patch("atomization_scorer.pipeline.true_pipeline.diagnose_paf_overlaps")
@patch("atomization_scorer.pipeline.true_pipeline.filter_paf")
@patch("atomization_scorer.pipeline.true_pipeline.align_with_minimap2")
@patch("atomization_scorer.pipeline.true_pipeline.extract_representatives")
def test_compute_true_alignment_propagates_filter_failure(
    mock_extract,
    mock_minimap2,
    mock_filter_paf,
    mock_diagnose,
    mock_resolve_paf,
    mock_validate_geese,
    mock_paf_to_geese,
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
    mock_diagnose.assert_not_called()
    mock_resolve_paf.assert_not_called()
    mock_paf_to_geese.assert_not_called()
    mock_validate_geese.assert_not_called()


# --------------------------------------------------------------------------------------
# Test: true pipeline returns after GEESE conversion
# --------------------------------------------------------------------------------------
@patch("atomization_scorer.pipeline.true_pipeline.paf_to_geese")
@patch("atomization_scorer.pipeline.true_pipeline.validate_non_overlapping_geese")
@patch("atomization_scorer.pipeline.true_pipeline.resolve_paf_overlaps")
@patch("atomization_scorer.pipeline.true_pipeline.diagnose_paf_overlaps")
@patch("atomization_scorer.pipeline.true_pipeline.filter_paf")
@patch("atomization_scorer.pipeline.true_pipeline.align_with_minimap2")
@patch("atomization_scorer.pipeline.true_pipeline.extract_representatives")
def test_compute_true_alignment_finishes_after_geese_conversion(
    mock_extract,
    mock_minimap2,
    mock_filter_paf,
    mock_diagnose,
    mock_resolve_paf,
    mock_validate_geese,
    mock_paf_to_geese,
    mini_fasta: Path,
    mini_geese: Path,
    output_dir: Path
):
    """compute_true_alignment should stop after generating the true GEESE file."""
    geese_path = compute_true_alignment(
        genomes_file=mini_fasta,
        atomization_file=mini_geese,
        output_directory=output_dir,
        run_overlap_diagnostics=False,
    )

    mock_extract.assert_called_once()
    mock_minimap2.assert_called_once()
    mock_filter_paf.assert_called_once()
    mock_diagnose.assert_not_called()
    mock_resolve_paf.assert_called_once()
    mock_paf_to_geese.assert_called_once()
    mock_validate_geese.assert_called_once()
    assert geese_path == output_dir / "true_atomization.geese"


# --------------------------------------------------------------------------------------
# Test: optionally runs overlap diagnostics between filter and resolve stages
# --------------------------------------------------------------------------------------
@patch("atomization_scorer.pipeline.true_pipeline.paf_to_geese")
@patch("atomization_scorer.pipeline.true_pipeline.validate_non_overlapping_geese")
@patch("atomization_scorer.pipeline.true_pipeline.resolve_paf_overlaps")
@patch("atomization_scorer.pipeline.true_pipeline.diagnose_paf_overlaps")
@patch("atomization_scorer.pipeline.true_pipeline.filter_paf")
@patch("atomization_scorer.pipeline.true_pipeline.align_with_minimap2")
@patch("atomization_scorer.pipeline.true_pipeline.extract_representatives")
def test_compute_true_alignment_optionally_runs_overlap_diagnostics(
    _mock_extract,
    _mock_minimap2,
    _mock_filter_paf,
    mock_diagnose,
    mock_resolve_paf,
    _mock_validate_geese,
    _mock_paf_to_geese,
    mini_fasta: Path,
    mini_geese: Path,
    output_dir: Path,
):
    """compute_true_alignment should optionally generate diagnostics from the filtered PAF."""
    compute_true_alignment(
        genomes_file=mini_fasta,
        atomization_file=mini_geese,
        output_directory=output_dir,
        run_overlap_diagnostics=True,
        overlap_report_min_len=200,
        overlap_plot_min_len=900,
    )

    mock_diagnose.assert_called_once_with(
        paf_file=output_dir / "minimap2_alignment_filtered.paf",
        representatives_fasta=output_dir / "mash_representatives.fa",
        output_directory=output_dir / "overlap_diagnostics",
        minimum_report_overlap_length=200,
        minimum_plot_overlap_length=900,
    )
    mock_resolve_paf.assert_called_once()
