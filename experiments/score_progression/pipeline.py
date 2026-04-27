"""
Run the full score progression pipeline: analyse -> extract_scores -> plot_scores.
"""

# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------
from __future__ import annotations

import logging
from pathlib import Path

from analysis import analyse
from data import extract_scores
from plots import plot_scores

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Score Progression Pipeline
# ---------------------------------------------------------------------------
def run_score_progression(
    fasta: Path,
    geese_bin: str,
    output_dir: Path,
    lastdb_threads: int = 16,
    lastal_threads: int = 15,
    min_length: int = 1000,
    min_alignment_length: int = 500,
    min_identity: int = 95,
    start_n: int = analyse.START_N,
    end_n: int | None = analyse.END_N,
    step: int = analyse.STEP,
    resume: bool = False,
) -> None:
    """Score each genome subset, extract score columns into a TSV, and generate a plot."""
    log.info("=" * 60)
    log.info("Step 1/3  analyse")
    results_csv = analyse.run(
        fasta=fasta,
        geese_bin=geese_bin,
        output_dir=output_dir,
        lastdb_threads=lastdb_threads,
        lastal_threads=lastal_threads,
        min_length=min_length,
        min_alignment_length=min_alignment_length,
        min_identity=min_identity,
        start_n=start_n,
        end_n=end_n,
        step=step,
        resume=resume,
    )

    log.info("=" * 60)
    log.info("Step 2/3  extract_scores")
    scores_tsv = extract_scores.run(
        input_csv=results_csv,
        output_tsv=output_dir / "scores.tsv",
    )

    log.info("=" * 60)
    log.info("Step 3/3  plot_scores")
    plot_scores.plot_scores(
        scores_tsv=scores_tsv,
        output_path=output_dir / "score_progression.png",
    )

    log.info("=" * 60)
    log.info("Pipeline complete.  outputs: %s", output_dir)


# ---------------------------------------------------------------------------
# Script Entry Point
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-8s  %(message)s",
        datefmt="%H:%M:%S",
    )
    args = analyse.parse_args()
    run_score_progression(
        fasta=args.fasta,
        geese_bin=args.geese,
        output_dir=args.output,
        lastdb_threads=args.lastdb_threads,
        lastal_threads=args.lastal_threads,
        min_length=args.min_length,
        min_alignment_length=args.min_aln_length,
        min_identity=args.min_ident,
        start_n=args.start_n,
        end_n=args.end_n,
        step=args.step,
        resume=args.resume,
    )
