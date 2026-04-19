"""
Read scores.tsv and generate a score-vs-genome-count line plot as a PNG file.
"""

# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------
from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
_HERE = Path(__file__).parent


# ---------------------------------------------------------------------------
# Command Line Interface
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Plot score progression from scores.tsv.")
    parser.add_argument(
        "--input",
        type=Path,
        default=_HERE / "results_genomes2_to582_step20" / "scores.tsv",
        help="Path to the scores TSV produced by extract_scores",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=_HERE / "results_genomes2_to582_step20" / "score_progression.png",
        help="Path where the PNG plot will be saved",
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def plot_scores(scores_tsv: Path, output_path: Path) -> None:
    """Read scores.tsv and save a score progression plot as a PNG file."""
    if not scores_tsv.is_file():
        raise FileNotFoundError(f"Scores TSV not found: {scores_tsv}")

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        log.warning("matplotlib not installed -- skipping plot generation")
        return

    with open(scores_tsv, newline="") as file:
        rows = list(csv.DictReader(file, delimiter="\t"))

    n_genomes = []
    overall = []
    alignment = []
    coverage = []

    for row in rows:
        n_genomes.append(int(row["n_genomes"].strip()))
        overall.append(float(row["overall_score"].strip()))
        alignment.append(float(row["alignment_score"].strip()))
        coverage.append(float(row["coverage_score"].strip()))

    if not n_genomes:
        log.warning("No rows in scores TSV - skipping plot generation")
        return

    fig, ax = plt.subplots(figsize=(max(10.0, len(n_genomes) * 0.8), 6.0))

    ax.plot(n_genomes, overall, "o-", label="Overall", linewidth=2, markersize=6, color="#2196F3")
    ax.plot(n_genomes, alignment, "s--", label="Alignment", linewidth=1.5, markersize=5, color="#FF9800")
    ax.plot(n_genomes, coverage, "^--", label="Coverage", linewidth=1.5, markersize=5, color="#4CAF50")

    ax.set_xlabel("Number of genomes")
    ax.set_ylabel("Score")
    ax.set_title("Score Progression as Genomes Are Added")
    ax.set_ylim(0, 1.05)
    ax.set_xticks(n_genomes)
    ax.legend(loc="lower left")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    log.info("Plot saved to %s", output_path)


# ---------------------------------------------------------------------------
# Script Entry Point
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    plot_scores(args.input, args.output)
