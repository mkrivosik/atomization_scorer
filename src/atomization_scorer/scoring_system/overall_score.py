"""
overall_score.py

Provides the main overall scoring function for genome atomization.

Modules
-------
compute_overall_score   : Computes overall atomization score.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
from __future__ import annotations
import logging
from pathlib import Path

from atomization_scorer.scoring_system import compute_alignment_score, compute_coverage_score

# --------------------------------------------------------------------------------------
# Logging
# --------------------------------------------------------------------------------------
logger = logging.getLogger(__name__)

# --------------------------------------------------------------------------------------
# Overall Score Function
# --------------------------------------------------------------------------------------

def compute_overall_score(
    genomes_file: Path,
    atomization_file: Path,
    output_directory: Path,
    level: str = "interval",
    per_class: bool = False,
    min_overlap_ratio: float = 0.8,
    alignment_weight: float = 0.7,
    coverage_weight: float = 0.3
) -> float:
    """
    Computes the overall atomization score using alignment and coverage scores.

    Parameters
    ----------
    genomes_file : Path
        Input genomes FASTA file containing the genome sequences.
    atomization_file : Path
        Input GEESE file containing the predicted atomization.
    output_directory : Path
        Path to the output directory where results are stored.
    level : str, optional, default: "interval"
        Alignment-score evaluation level.
    per_class : bool, optional, default: False
        Whether to compute per-class alignment metrics.
    min_overlap_ratio : float, optional, default: 0.8
        Minimum overlap ratio for interval-level alignment scoring.
    alignment_weight : float, optional, default: 0.7
        Weight of the alignment score in the geometric mean.
    coverage_weight : float, optional, default: 0.3
        Weight of the coverage score in the geometric mean.

    Raises
    ------
    FileNotFoundError
        Raised if genomes_file or atomization_file do not exist.
    ValueError
        Raised if the weights are negative or do not sum to 1.0.

    Returns
    -------
    float
        The overall atomization score.
    """
    if not genomes_file.is_file():
        raise FileNotFoundError(f"Genomes FASTA file not found: {genomes_file}")
    if not atomization_file.is_file():
        raise FileNotFoundError(f"Atomization file not found: {atomization_file}")
    if alignment_weight < 0 or coverage_weight < 0:
        raise ValueError("Alignment and coverage weights must be non-negative.")
    if alignment_weight + coverage_weight != 1.0:
        raise ValueError("Alignment and coverage weights must sum to 1.0.")

    output_directory.mkdir(parents=True, exist_ok=True)

    logger.info(
        "Computing overall score for genomes=%s, atomization=%s, level=%s, per_class=%s, "
        "min_overlap_ratio=%s, alignment_weight=%s, coverage_weight=%s",
        genomes_file,
        atomization_file,
        level,
        per_class,
        min_overlap_ratio,
        alignment_weight,
        coverage_weight,
    )

    logger.info("Computing alignment score")
    alignment_score = compute_alignment_score(
        genomes_file=genomes_file,
        atomization_file=atomization_file,
        output_directory=output_directory,
        level=level,
        per_class=per_class,
        min_overlap_ratio=min_overlap_ratio
    )
    logger.info("Alignment score: %s", alignment_score)

    logger.info("Computing coverage score")
    coverage_score = compute_coverage_score(
        genomes_file=genomes_file,
        atomization_file=atomization_file
    )
    logger.info("Coverage score: %s", coverage_score)

    # Weighted geometric mean
    overall_score = (alignment_score ** alignment_weight) * (coverage_score ** coverage_weight)
    overall_score = min(max(overall_score, 0.0), 1.0)

    logger.info("Overall score: %s", overall_score)
    return overall_score
