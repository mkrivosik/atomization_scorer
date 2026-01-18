"""
overall_score.py

Provides the main overall scoring function for genome atomization.

Modules
-------
compute_overall_score : Computes overall atomization score.
"""

# ---------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------
from pathlib import Path
from atomization_scorer.scoring_system import compute_alignment_score
from atomization_scorer.scoring_system import compute_coverage_score

# ---------------------------------------------------------------------
# Overall Score Function
# ---------------------------------------------------------------------

def compute_overall_score(
    genomes_file: Path,
    atomization_file: Path,
    output_directory: Path
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

    Raises
    ------
    FileNotFoundError
        Raised if genomes_file or atomization_file do not exist.

    Returns
    -------
    float
        The overall atomization score.
    """
    if not genomes_file.is_file():
        raise FileNotFoundError(f"Genomes FASTA file not found: {genomes_file}")
    if not atomization_file.is_file():
        raise FileNotFoundError(f"Atomization file not found: {atomization_file}")

    output_directory.mkdir(parents=True, exist_ok=True)

    print("Computing alignment score...")
    alignment_score = compute_alignment_score(
        genomes_file=genomes_file,
        atomization_file=atomization_file,
        output_directory=output_directory,
        level="interval",
        per_class=False,
        min_overlap_ratio=0.8
    )
    print("Alignment score:", alignment_score)

    print("Computing coverage score...")
    coverage_score = compute_coverage_score(
        genomes_file=genomes_file,
        atomization_file=atomization_file
    )
    print("Coverage score:", coverage_score)

    # Weighted geometric mean
    overall_score = (alignment_score ** 0.7) * (coverage_score ** 0.3)
    return min(max(overall_score, 0.0), 1.0)
