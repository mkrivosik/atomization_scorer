"""
true_pipeline.py

Provides the pipeline used to compute true (gold standard) genome atomization.

Functions
---------
compute_true_alignment  : Compute true (gold standard) atomization and return a GEESE file path.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
from __future__ import annotations
import logging
from pathlib import Path

from atomization_scorer.data_processing import (
    filter_paf,
    paf_to_geese,
    resolve_paf_overlaps,
    validate_non_overlapping_geese,
)
from atomization_scorer.diagnostics import diagnose_paf_overlaps

from .minimap2_aligner import align_with_minimap2
from .representatives_selector import extract_representatives

# --------------------------------------------------------------------------------------
# Logging
# --------------------------------------------------------------------------------------
log = logging.getLogger(__name__)


# --------------------------------------------------------------------------------------
# True (Gold Standard) Alignment Pipeline
# --------------------------------------------------------------------------------------
def compute_true_alignment(
    genomes_file: Path,
    atomization_file: Path,
    output_directory: Path,
    mode: str = "mash",
    minimap2_preset: str = "asm20",
    minimap2_secondary_ratio: float = 0.1,
    minimap2_emit_cigar: bool = True,
    minimum_similarity: float = 0.95,
    minimum_alignment_length: int = 500,
    run_overlap_diagnostics: bool = False,
    diagnostics_directory: Path | None = None,
    minimum_report_overlap_length: int = 0,
    minimum_plot_overlap_length: int = 0,
    overlap_include_reverse: bool = False,
    run_dotter: bool = True,
    dotter_extra_args: list[str] | None = None,
) -> Path:
    """
    Run a full true (gold standard) genome atomization pipeline: extract representatives, align
    genome sequences on them with minimap2, filter PAF alignments, resolve PAF overlaps, convert
    to GEESE format, and validate the result.

    Parameters
    ----------
    genomes_file : Path
        Input genomes FASTA file containing the genomes sequences.
    atomization_file : Path
        Input GEESE file containing the predicted atomization.
    output_directory : Path
        Path to the output directory where results are stored.
    mode : str, optional, default="mash"
        Representative selection mode ("mash" or "first").
    minimap2_preset : str, optional, default="asm20"
        Minimap2 preset passed to -x.
    minimap2_secondary_ratio : float, optional, default=0.1
        Minimap2 secondary-to-primary score ratio passed to -p.
    minimap2_emit_cigar : bool, optional, default=True
        Whether to include CIGAR strings in PAF output (-c flag).
    minimum_similarity : float, optional, default=0.95
        Minimum similarity for PAF filtering.
    minimum_alignment_length : int, optional, default=500
        Minimum alignment length for PAF filtering.
    run_overlap_diagnostics : bool, optional, default=False
        Whether to generate overlap-diagnostic reports from the filtered PAF before
        overlap resolution.
    diagnostics_directory : Path or None, optional, default=None
        Directory where overlap diagnostics should be written. Defaults to
        output_directory / "overlap_diagnostics" when not provided.
    minimum_report_overlap_length : int, optional, default=0
        Minimum overlap length required for overlap-level reporting.
    minimum_plot_overlap_length : int, optional, default=0
        Minimum overlap length required for dotplot FASTA generation.
    overlap_include_reverse : bool, optional, default=False
        Whether to duplicate anchor-partner diagnostics from both perspectives.
    run_dotter : bool, optional, default=True
        Whether to run Dotter immediately after generating anchor FASTA inputs.
    dotter_extra_args : list[str] or None, optional, default=None
        Additional command-line arguments passed to Dotter when run_dotter is enabled.

    Raises
    ------
    FileNotFoundError
        Raised if genomes_file or atomization_file do not exist.

    Returns
    -------
    Path
        Generated true (gold standard) GEESE file path.
    """
    if not genomes_file.is_file():
        raise FileNotFoundError(f"Genomes FASTA file not found: {genomes_file}")
    if not atomization_file.is_file():
        raise FileNotFoundError(f"Atomization file not found: {atomization_file}")

    output_directory.mkdir(parents=True, exist_ok=True)
    log.info("=" * 60)
    log.info(
        "Computing true alignment with mode=%s for genomes=%s atomization=%s into output=%s",
        mode,
        genomes_file,
        atomization_file,
        output_directory,
    )

    # Extract representatives
    log.info("=" * 60)
    representatives_fasta = output_directory / f"{mode}_representatives.fa"
    log.info("Extracting representatives to %s", representatives_fasta)
    extract_representatives(
        genomes_file=genomes_file,
        atomization_file=atomization_file,
        output_path=representatives_fasta,
        mode=mode
    )

    # Minimap2 alignment
    log.info("=" * 60)
    paf_file = output_directory / "minimap2_alignments.paf"
    log.info("Running minimap2 alignment to %s", paf_file)
    align_with_minimap2(
        query=genomes_file,
        target=representatives_fasta,
        output_path=paf_file,
        preset=minimap2_preset,
        secondary_ratio=minimap2_secondary_ratio,
        emit_cigar=minimap2_emit_cigar,
    )

    # Filter PAF
    log.info("=" * 60)
    filtered_paf = output_directory / "minimap2_alignment_filtered.paf"
    log.info(
        "Filtering PAF into %s with minimum_similarity=%s and minimum_alignment_length=%s",
        filtered_paf,
        minimum_similarity,
        minimum_alignment_length,
    )
    filter_paf(
        paf_file=paf_file,
        output_file=filtered_paf,
        minimum_similarity=minimum_similarity,
        minimum_alignment_length=minimum_alignment_length
    )

    # Run Diagnostics
    if run_overlap_diagnostics:
        log.info("=" * 60)
        resolved_diagnostics_directory = (
            diagnostics_directory
            if diagnostics_directory is not None
            else output_directory / "overlap_diagnostics"
        )
        log.info(
            "Generating overlap diagnostics into %s with minimum_report_overlap_length=%s and minimum_plot_overlap_length=%s",
            resolved_diagnostics_directory,
            minimum_report_overlap_length,
            minimum_plot_overlap_length,
        )
        diagnose_paf_overlaps(
            paf_file=filtered_paf,
            representatives_fasta=representatives_fasta,
            output_directory=resolved_diagnostics_directory,
            minimum_report_overlap_length=minimum_report_overlap_length,
            minimum_plot_overlap_length=minimum_plot_overlap_length,
            include_reverse=overlap_include_reverse,
            run_dotter=run_dotter,
            dotter_extra_args=dotter_extra_args,
        )

    # Resolve PAF overlaps
    log.info("=" * 60)
    resolved_paf = output_directory / "minimap2_alignment_resolved.paf"
    log.info("Resolving overlapping PAF alignments into %s", resolved_paf)
    resolve_paf_overlaps(
        paf_file=filtered_paf,
        output_file=resolved_paf,
        minimum_alignment_length=minimum_alignment_length,
    )

    # Convert PAF to GEESE
    log.info("=" * 60)
    geese_file = output_directory / "true_atomization.geese"
    log.info("Converting resolved PAF to GEESE at %s", geese_file)
    paf_to_geese(
        paf_file=resolved_paf,
        output_file=geese_file
    )

    # Validate true atomization
    log.info("=" * 60)
    log.info("Validating true atomization non-overlap at %s", geese_file)
    validate_non_overlapping_geese(geese_file)

    log.info("=" * 60)
    log.info("True alignment pipeline finished with output %s", geese_file)
    return geese_file
