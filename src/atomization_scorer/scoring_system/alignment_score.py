"""
alignment_score.py

Computes the alignment-based atomization score.

Functions
---------
compute_alignment_score : Computes the alignment score comparing predicted
                          atomization to gold-standard (true) atomization.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
from __future__ import annotations
import logging
from pathlib import Path

from atomization_scorer.pipeline import compute_true_alignment
from atomization_scorer.visualization import plot_atomization

from .base_metrics import compute_base_level_metrics
from .interval_metrics import compute_interval_level_metrics

# --------------------------------------------------------------------------------------
# Logging
# --------------------------------------------------------------------------------------
log = logging.getLogger(__name__)


# --------------------------------------------------------------------------------------
# Alignment Score Function
# --------------------------------------------------------------------------------------
def compute_alignment_score(
    genomes_file: Path,
    atomization_file: Path,
    output_directory: Path,
    level: str = "interval",
    per_class: bool = False,
    minimum_overlap_ratio: float = 0.8,
    representative_mode: str = "mash",
    minimum_similarity: float = 0.95,
    minimum_alignment_length: int = 500,
    minimap2_preset: str = "asm20",
    minimap2_secondary_ratio: float = 0.1,
    minimap2_emit_cigar: bool = True,
    run_overlap_diagnostics: bool = False,
    minimum_report_overlap_length: int = 0,
    minimum_plot_overlap_length: int = 0,
    overlap_include_reverse: bool = False,
    run_dotter: bool = True,
) -> float | list[dict[str, int | float]]:
    """
    Compute alignment score (F1-score) at base or interval level, optionally per class.

    Parameters
    ----------
    genomes_file : Path
        Input genomes FASTA file containing the genome sequences.
    atomization_file : Path
        Input GEESE file containing the predicted atomization.
    output_directory : Path
        Path to the output directory where results are stored.
    level : str, optional, default="interval"
        Select "base" or "interval".
    per_class : bool, optional, default=False
        If True, compute score per class; else overall.
    minimum_overlap_ratio : float, optional, default=0.8
        Minimum overlap ratio for interval-level scoring.
    representative_mode : str, optional, default="mash"
        Representative selection strategy: "mash" or "first".
    minimum_similarity : float, optional, default=0.95
        Minimum similarity threshold for PAF filtering.
    minimum_alignment_length : int, optional, default=500
        Minimum alignment length threshold for PAF filtering.
    minimap2_preset : str, optional, default="asm20"
        Minimap2 preset passed to -x.
    minimap2_secondary_ratio : float, optional, default=0.1
        Minimap2 secondary-to-primary score ratio passed to -p.
    minimap2_emit_cigar : bool, optional, default=True
        Whether to include CIGAR strings in PAF output (-c flag).
    run_overlap_diagnostics : bool, optional, default=False
        Whether to generate overlap-diagnostic reports before overlap resolution.
    minimum_report_overlap_length : int, optional, default=0
        Minimum overlap length for overlap-level reporting.
    minimum_plot_overlap_length : int, optional, default=0
        Minimum overlap length for dotplot FASTA generation.
    overlap_include_reverse : bool, optional, default=False
        Whether to duplicate anchor-partner diagnostics from both perspectives.
    run_dotter : bool, optional, default=True
        Whether to run Dotter immediately after generating anchor FASTA inputs.

    Raises
    ------
    FileNotFoundError
        Raised if genomes_file or atomization_file do not exist.
    ValueError
        If the level is not "base" or "interval".

    Returns
    -------
    float or List[Dict[str, int | float]]
        If per_class is False, returns overall interval-level F1-score between 0.0 and 1.0.
        If per_class is True, returns a list of dictionaries, each containing:
            "Class": int -> atomization class,
            "F1-score": float -> F1-score for that class (0.0 to 1.0).
    """
    if not genomes_file.is_file():
        raise FileNotFoundError(f"Genomes FASTA file not found: {genomes_file}")
    if not atomization_file.is_file():
        raise FileNotFoundError(f"Atomization file not found: {atomization_file}")

    if level not in ("base", "interval"):
        raise ValueError("Level must be 'base' or 'interval'.")

    output_directory.mkdir(parents=True, exist_ok=True)
    log.info(
        "Computing alignment score with level=%s per_class=%s minimum_overlap_ratio=%s",
        level,
        per_class,
        minimum_overlap_ratio,
    )

    pipeline_directory = output_directory / "true_atomization"
    diagnostics_directory = output_directory / "overlap_diagnostics"
    metrics_directory = output_directory / "metrics"
    visualization_directory = output_directory / "visualization"

    log.info("=" * 60)
    log.info("Computing gold standard (true) alignment")
    true_geese = compute_true_alignment(
        genomes_file=genomes_file,
        atomization_file=atomization_file,
        output_directory=pipeline_directory,
        mode=representative_mode,
        minimum_similarity=minimum_similarity,
        minimum_alignment_length=minimum_alignment_length,
        minimap2_preset=minimap2_preset,
        minimap2_secondary_ratio=minimap2_secondary_ratio,
        minimap2_emit_cigar=minimap2_emit_cigar,
        run_overlap_diagnostics=run_overlap_diagnostics,
        minimum_report_overlap_length=minimum_report_overlap_length,
        minimum_plot_overlap_length=minimum_plot_overlap_length,
        overlap_include_reverse=overlap_include_reverse,
        run_dotter=run_dotter,
        diagnostics_directory=diagnostics_directory,
    )

    log.info("=" * 60)
    if level == "base":
        log.info("Computing base-level metrics")
        score = compute_base_level_metrics(
            predicted_geese=atomization_file,
            true_geese=Path(true_geese),
            output_directory=metrics_directory,
            per_class=per_class
        )
    else:
        log.info("Computing interval-level metrics")
        score = compute_interval_level_metrics(
            predicted_geese=atomization_file,
            true_geese=Path(true_geese),
            output_directory=metrics_directory,
            per_class=per_class,
            minimum_overlap_ratio=minimum_overlap_ratio
        )

    log.info("=" * 60)
    log.info("Generating atomization visualization into %s", visualization_directory)
    plot_atomization(
        genomes_file=genomes_file,
        true_atoms_file=Path(true_geese),
        predicted_atoms_file=atomization_file,
        output_directory=visualization_directory,
    )

    log.info("=" * 60)
    log.info("Alignment score result: %s", score)
    return score
