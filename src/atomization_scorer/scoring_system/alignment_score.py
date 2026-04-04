"""
alignment_score.py

Computes the alignment-based atomization score.

Modules
-------
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
logger = logging.getLogger(__name__)

# --------------------------------------------------------------------------------------
# Alignment Score Function
# --------------------------------------------------------------------------------------

def compute_alignment_score(
    genomes_file: Path,
    atomization_file: Path,
    output_directory: Path,
    level: str = "interval",
    per_class: bool = False,
    min_overlap_ratio: float = 0.8
) -> float | list[dict[str, int | float]]:
    """
    Compute alignment score (F1-score) at base or interval level, optionally per genome.

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
    min_overlap_ratio : float, optional, default=0.8
        Minimum overlap ratio for interval-level scoring.

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
    logger.info(
        "Computing alignment score with level=%s per_class=%s min_overlap_ratio=%s",
        level,
        per_class,
        min_overlap_ratio,
    )

    logger.info("Computing gold standard (true) alignment")
    true_geese = compute_true_alignment(
        genomes_file=genomes_file,
        atomization_file=atomization_file,
        output_directory=output_directory
    )

    visualization_directory = output_directory / "atomization_visualizations"
    if level == "base":
        logger.info("Computing base-level metrics")
        score = compute_base_level_metrics(
            predicted_geese=atomization_file,
            true_geese=Path(true_geese),
            output_directory=output_directory,
            per_class=per_class
        )
    else:
        logger.info("Computing interval-level metrics")
        score = compute_interval_level_metrics(
            predicted_geese=atomization_file,
            true_geese=Path(true_geese),
            output_directory=output_directory,
            per_class=per_class,
            min_overlap_ratio=min_overlap_ratio
        )

    logger.info("Generating atomization visualization into %s", visualization_directory)
    plot_atomization(
        genomes_file=genomes_file,
        true_atoms_file=Path(true_geese),
        predicted_atoms_file=atomization_file,
        output_directory=visualization_directory,
    )

    logger.info("Alignment score result: %s", score)
    return score
