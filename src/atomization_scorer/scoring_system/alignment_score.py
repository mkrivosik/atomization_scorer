"""
alignment_score.py

Computes the alignment-based atomization score.

Modules
-------
compute_alignment_score : Computes the alignment score comparing predicted
                          atomization to gold-standard (true) atomization.
"""

# ---------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------
from pathlib import Path

from atomization_scorer.pipeline import compute_true_alignment

from .base_metrics import compute_base_level_metrics
from .interval_metrics import compute_interval_level_metrics

# ---------------------------------------------------------------------
# Alignment Score Function
# ---------------------------------------------------------------------

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
        Input genomes FASTA file containing the genomes sequences.
    atomization_file : Path
        Input GEESE file containing the predicted atomization.
    output_directory : Path
        Path to the output directory where results are stored.
    level : str, optional, default: "interval"
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
        If level is not "base" or "interval".

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

    output_directory.mkdir(parents=True, exist_ok=True)

    print("Computing gold standard (true) alignment...")
    true_geese = compute_true_alignment(genomes_file, atomization_file, output_directory)

    print(
        f"Computing score with parameters\n"
        f"  level: {level}\n"
        f"  per_class: {per_class}\n"
        f"  min_overlap_ratio: {min_overlap_ratio}\n"
    )
    if level == "base":
        score = compute_base_level_metrics(
            predicted_geese=atomization_file,
            true_geese=true_geese,
            output_directory=output_directory,
            per_class=per_class
        )
    elif level == "interval":
        score = compute_interval_level_metrics(
            predicted_geese=atomization_file,
            true_geese=true_geese,
            output_directory=output_directory,
            per_class=per_class,
            min_overlap_ratio=min_overlap_ratio
        )
    else:
        raise ValueError("Level must be 'base' or 'interval'.")

    return score
