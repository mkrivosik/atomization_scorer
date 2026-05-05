"""
overall_score.py

Provides the main overall scoring function for genome atomization.

Functions
---------
compute_overall_score   : Computes overall atomization score.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
from __future__ import annotations
import logging
import math
from pathlib import Path

from atomization_scorer.scoring_system import compute_alignment_score, compute_coverage_score

# --------------------------------------------------------------------------------------
# Logging
# --------------------------------------------------------------------------------------
log = logging.getLogger(__name__)


# --------------------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------------------
def _write_scores(
    output_directory: Path,
    alignment_score: float,
    coverage_score: float,
    overall_score: float,
    alignment_weight: float,
    coverage_weight: float,
) -> None:
    scores_file = output_directory / "scores.tsv"
    width_metric = 25
    rows = [
        ("alignment_score",  f"{alignment_score:.6f}"),
        ("coverage_score",   f"{coverage_score:.6f}"),
        ("overall_score",    f"{overall_score:.6f}"),
        ("alignment_weight", f"{alignment_weight:.6f}"),
        ("coverage_weight",  f"{coverage_weight:.6f}"),
    ]
    with scores_file.open("w") as file:
        file.write(f"{'metric':<{width_metric}} value\n")
        file.write(f"{'-' * width_metric} {'-' * 8}\n")
        for metric, value in rows:
            file.write(f"{metric:<{width_metric}} {value}\n")


def _write_scores_per_class(
    output_directory: Path,
    results: list[dict[str, int | float]],
    alignment_weight: float,
    coverage_weight: float,
) -> None:
    scores_file = output_directory / "scores_per_class.tsv"
    width_class, width_alignment_score, width_coverage_score, width_overall = 8, 16, 16, 12
    width_metric = 25
    with scores_file.open("w") as file:
        file.write(
            f"{'class':<{width_class}} "
            f"{'alignment_score':<{width_alignment_score}} "
            f"{'coverage_score':<{width_coverage_score}} "
            f"{'overall'}\n"
        )
        file.write(
            f"{'-' * width_class} "
            f"{'-' * width_alignment_score} "
            f"{'-' * width_coverage_score} "
            f"{'-' * width_overall}\n"
        )
        for class_entry in results:
            file.write(
                f"{class_entry['Class']:<{width_class}} "
                f"{class_entry['Alignment_score']:<{width_alignment_score}.6f} "
                f"{class_entry['Coverage_score']:<{width_coverage_score}.6f} "
                f"{class_entry['Overall']:.6f}\n"
            )
        file.write("\n")
        file.write(f"{'alignment_weight':<{width_metric}} {alignment_weight:.6f}\n")
        file.write(f"{'coverage_weight':<{width_metric}} {coverage_weight:.6f}\n")


# --------------------------------------------------------------------------------------
# Overall Score Function
# --------------------------------------------------------------------------------------
def compute_overall_score(
    genomes_file: Path,
    atomization_file: Path,
    output_directory: Path,
    level: str = "interval",
    per_class: bool = False,
    minimum_overlap_ratio: float = 0.8,
    alignment_weight: float = 0.7,
    coverage_weight: float = 0.3,
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
    Computes the overall atomization score using alignment and coverage scores.

    Parameters
    ----------
    genomes_file : Path
        Input genomes FASTA file containing the genome sequences.
    atomization_file : Path
        Input GEESE file containing the predicted atomization.
    output_directory : Path
        Path to the output directory where results are stored.
    level : str, optional, default="interval"
        Alignment score evaluation level ("base" or "interval").
    per_class : bool, optional, default=False
        Whether to compute per-class alignment and coverage metrics and return per-class scores.
    minimum_overlap_ratio : float, optional, default=0.8
        Minimum overlap ratio for interval-level alignment scoring.
    alignment_weight : float, optional, default=0.7
        Weight of the alignment score in the geometric mean.
    coverage_weight : float, optional, default=0.3
        Weight of the coverage score in the geometric mean.
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
        Raised if the weights are negative or do not sum to 1.0.

    Returns
    -------
    float or list[dict[str, int | float]]
        If per_class is False, returns the overall atomization score (0.0 to 1.0).
        If per_class is True, returns a list of dictionaries sorted by class, each containing:
            "Class": int -> atomization class,
            "Alignment_score": float -> alignment F1-score for that class,
            "Coverage_score": float -> coverage fraction for that class,
            "Overall": float -> weighted geometric mean for that class.
    """
    if not genomes_file.is_file():
        raise FileNotFoundError(f"Genomes FASTA file not found: {genomes_file}")
    if not atomization_file.is_file():
        raise FileNotFoundError(f"Atomization file not found: {atomization_file}")
    if alignment_weight < 0 or coverage_weight < 0:
        raise ValueError("Alignment and coverage weights must be non-negative.")
    if not math.isclose(alignment_weight + coverage_weight, 1.0):
        raise ValueError("Alignment and coverage weights must sum to 1.0.")

    output_directory.mkdir(parents=True, exist_ok=True)

    log.info(
        "Computing overall score for genomes=%s, atomization=%s, level=%s, per_class=%s, "
        "minimum_overlap_ratio=%s, alignment_weight=%s, coverage_weight=%s",
        genomes_file,
        atomization_file,
        level,
        per_class,
        minimum_overlap_ratio,
        alignment_weight,
        coverage_weight,
    )

    log.info("=" * 60)
    log.info("Computing alignment score")
    alignment_score = compute_alignment_score(
        genomes_file=genomes_file,
        atomization_file=atomization_file,
        output_directory=output_directory,
        level=level,
        per_class=per_class,
        minimum_overlap_ratio=minimum_overlap_ratio,
        representative_mode=representative_mode,
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
    )
    log.info("Alignment score: %s", alignment_score)

    log.info("=" * 60)
    log.info("Computing coverage score")
    coverage_score = compute_coverage_score(
        genomes_file=genomes_file,
        atomization_file=atomization_file,
        per_class=per_class,
    )
    log.info("Coverage score: %s", coverage_score)

    log.info("=" * 60)
    log.info("Computing overall score")

    if per_class:
        assert isinstance(alignment_score, list)
        assert isinstance(coverage_score, list)
        alignment_by_class = {entry["Class"]: entry["F1-score"] for entry in alignment_score}
        coverage_by_class = {entry["Class"]: entry["Coverage"] for entry in coverage_score}
        all_classes = sorted(set(alignment_by_class) | set(coverage_by_class))

        results = []
        for atom_class in all_classes:
            f1_score = alignment_by_class.get(atom_class, 0.0)
            coverage = coverage_by_class.get(atom_class, 0.0)
            overall = min(max((f1_score ** alignment_weight) * (coverage ** coverage_weight), 0.0), 1.0)
            results.append({"Class": atom_class, "Alignment_score": f1_score, "Coverage_score": coverage, "Overall": overall})

        _write_scores_per_class(output_directory, results, alignment_weight, coverage_weight)
        log.info("Per-class overall scores: %s", results)
        return results

    assert isinstance(alignment_score, float)
    assert isinstance(coverage_score, float)
    overall_score = min(max((alignment_score ** alignment_weight) * (coverage_score ** coverage_weight), 0.0), 1.0)

    _write_scores(
        output_directory=output_directory,
        alignment_score=alignment_score,
        coverage_score=coverage_score,
        overall_score=overall_score,
        alignment_weight=alignment_weight,
        coverage_weight=coverage_weight,
    )

    log.info("Overall score: %s", overall_score)
    return overall_score
