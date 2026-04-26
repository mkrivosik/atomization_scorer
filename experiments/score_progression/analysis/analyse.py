"""
Incrementally score genome subsets and record how the atomization score changes as genomes are added.
"""

# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------
from __future__ import annotations

import argparse
import csv
import logging
import subprocess
import sys
from pathlib import Path
from typing import cast

from atomization_scorer import (
    compute_alignment_score,
    compute_coverage_score,
    read_fasta,
    write_fasta,
)

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
_HERE        = Path(__file__).parent
FASTA        = _HERE.parent.parent.parent / "tests" / "fixtures" / "mini.fa"
GEESE_BIN    = Path.home() / "bin/geese"

START_N      = 2    # first genome count to score (must be >= 2)
END_N        = 596
STEP         = 20

OUTPUT_DIR   = _HERE.parent / f"results_genomes{START_N}_to{END_N}_step{STEP}"

CSV_FIELDS = [
    "n_genomes",
    "added_genome",
    "genome_list",
    "overall_score",
    "alignment_score",
    "coverage_score",
    "error",
]


# ---------------------------------------------------------------------------
# External Tool Runners
# ---------------------------------------------------------------------------
def run_lastdb(fasta: Path, db_prefix: Path, threads: int = 16) -> None:
    """Build a LAST sequence database from a FASTA file."""
    cmd = ["lastdb", "-v", "-P", str(threads), str(db_prefix), str(fasta)]
    log.info("  lastdb: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise subprocess.CalledProcessError(
            result.returncode, cmd, output=result.stdout, stderr=result.stderr
        )


def run_lastal(db_prefix: Path, fasta: Path, output_psl: Path, threads: int = 15) -> None:
    """Align a FASTA file against a LAST database and write the result as a PSL file."""
    cmd = (
        f"lastal -v -P {threads} {db_prefix} {fasta} "
        f"| maf-convert psl > {output_psl}"
    )
    log.info("  lastal: %s", cmd)
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise subprocess.CalledProcessError(
            result.returncode, cmd, output=result.stdout, stderr=result.stderr
        )


def run_geese(
    geese_bin: str,
    psl: Path,
    output_geese: Path,
    *,
    min_length: int = 1000,
    min_alignment_length: int = 500,
    min_identity: int = 95,
) -> None:
    """Run GEESE on a PSL file to produce a predicted atomization."""
    cmd = (
        f"{geese_bin} {psl} "
        f"--minLength {min_length} "
        f"--minAlnLength {min_alignment_length} "
        f"--minIdent {min_identity} "
        f"> {output_geese}"
    )
    log.info("  geese:  %s", cmd)
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise subprocess.CalledProcessError(
            result.returncode, cmd, output=result.stdout, stderr=result.stderr
        )

    if output_geese.stat().st_size == 0:
        raise RuntimeError(f"GEESE produced an empty file: {output_geese}")


# ---------------------------------------------------------------------------
# Atomization Scorer
# ---------------------------------------------------------------------------
def run_scorer(fasta: Path, geese: Path, output_dir: Path) -> dict[str, float | None]:
    """Score predicted atomization against the true (gold standard) atomization."""
    output_dir.mkdir(parents=True, exist_ok=True)

    alignment = None
    coverage = None

    try:
        alignment = cast(float, compute_alignment_score(
            genomes_file=fasta,
            atomization_file=geese,
            output_directory=output_dir,
            level="interval",
            per_class=False,
            min_overlap_ratio=0.8,
        ))
    except Exception as alignment_exc:
        log.exception("  Alignment score failed: %s", alignment_exc)

    try:
        coverage = compute_coverage_score(
            genomes_file=fasta,
            atomization_file=geese,
        )
    except Exception as coverage_exc:
        log.exception("  Coverage score failed: %s", coverage_exc)

    overall = None
    if alignment is not None and coverage is not None:
        overall = (alignment ** 0.7) * (coverage ** 0.3)
        overall = min(max(overall, 0.0), 1.0)

    return {"overall": overall, "alignment": alignment, "coverage": coverage}


# ---------------------------------------------------------------------------
# Results I/O
# ---------------------------------------------------------------------------
def write_results_csv(rows: list[dict], path: Path) -> None:
    """Write the accumulated results to a CSV file, overwriting any existing content."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=CSV_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def read_results_csv(path: Path) -> list[dict]:
    """Read an existing results CSV into a list of row dicts, or return [] if it does not exist."""
    if not path.is_file():
        return []
    with open(path, newline="") as file:
        return list(csv.DictReader(file))


# ---------------------------------------------------------------------------
# Command Line Interface
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Score Progression Analysis: incrementally add genomes "
            "and track how the atomization score changes."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "fasta",
        type=Path,
        nargs="?",
        default=FASTA,
        help="Input FASTA file containing all N genomes",
    )
    parser.add_argument(
        "--geese",
        default=str(GEESE_BIN),
        help="Path to the GEESE binary",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=OUTPUT_DIR,
        help="Output directory for per-iteration data and CSV",
    )
    parser.add_argument("--lastdb-threads", type=int, default=16)
    parser.add_argument("--lastal-threads", type=int, default=15)
    parser.add_argument("--min-length", type=int, default=1000, help="GEESE --minLength")
    parser.add_argument("--min-aln-length", type=int, default=500, help="GEESE --minAlnLength")
    parser.add_argument("--min-ident", type=int, default=95, help="GEESE --minIdent")
    parser.add_argument(
        "--start-n",
        type=int,
        default=START_N,
        help="Start from this many genomes",
    )
    parser.add_argument(
        "--end-n",
        type=int,
        default=END_N,
        help="Stop at this many genomes (inclusive)",
    )
    parser.add_argument(
        "--step",
        type=int,
        default=STEP,
        help="Increment genome count by this many each iteration",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Skip iterations whose results already exist in results.csv",
    )

    return parser.parse_args()


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
def run(
    fasta: Path,
    geese_bin: str,
    output_dir: Path,
    lastdb_threads: int = 16,
    lastal_threads: int = 15,
    min_length: int = 1000,
    min_alignment_length: int = 500,
    min_identity: int = 95,
    start_n: int = START_N,
    end_n: int | None = END_N,
    step: int = STEP,
    resume: bool = False,
) -> Path:
    """Iterate over genome subsets (start_n...end_n), score each, and write results.csv."""
    log.info("Reading genomes from %s", fasta)
    all_genomes = read_fasta(fasta)
    genome_names = list(all_genomes.keys())
    n_total = len(genome_names)
    if n_total == 0:
        raise ValueError(f"No sequences found in {fasta}")
    log.info("Found %d genomes: %s", n_total, ", ".join(genome_names))

    resolved_end_n = end_n if end_n is not None else n_total
    resolved_end_n = min(resolved_end_n, n_total)

    if start_n < 2:
        raise ValueError(f"start_n must be >= 2 (need at least 2 genomes), got {start_n}")
    if start_n > resolved_end_n:
        raise ValueError(f"start_n ({start_n}) > end_n ({resolved_end_n}), nothing to do")

    output_dir.mkdir(parents=True, exist_ok=True)
    results_csv = output_dir / "results_genomes2_to582_step20.csv"

    if resume:
        existing_results = read_results_csv(results_csv)
        completed_genome_counts = {int(result["n_genomes"]) for result in existing_results}
        results = list(existing_results)
        log.info("Resuming: %d iterations already completed", len(completed_genome_counts))
    else:
        completed_genome_counts = set()
        results = []

    steps = list(range(start_n, resolved_end_n + 1, step))
    if not steps or steps[-1] != resolved_end_n:
        steps.append(resolved_end_n)

    for n_genomes in steps:
        if n_genomes in completed_genome_counts:
            log.info("[n=%d] skipped (already completed)", n_genomes)
            continue

        subset_names = genome_names[:n_genomes]
        added_genome = genome_names[n_genomes - 1]
        iter_dir = output_dir / f"n{n_genomes:03d}_{added_genome}"
        iter_dir.mkdir(parents=True, exist_ok=True)

        log.info("=" * 60)
        log.info("[n=%d] Adding genome: %s", n_genomes, added_genome)
        log.info("[n=%d] Genomes: %s", n_genomes, ", ".join(subset_names))

        row = {
            "n_genomes": str(n_genomes),
            "added_genome": added_genome,
            "genome_list": ";".join(subset_names),
            "overall_score": "",
            "alignment_score": "",
            "coverage_score": "",
            "error": "",
        }

        try:
            # 1. Extract subset FASTA
            subset_fasta = iter_dir / "subset.fa"
            write_fasta({name: all_genomes[name] for name in subset_names}, subset_fasta)
            log.info("[n=%d] Wrote %d genomes to %s", n_genomes, n_genomes, subset_fasta)

            # 2. Run LAST alignment
            db_prefix = iter_dir / "atoms-db"
            run_lastdb(subset_fasta, db_prefix, threads=lastdb_threads)

            psl_file = iter_dir / "atoms.psl"
            run_lastal(db_prefix, subset_fasta, psl_file, threads=lastal_threads)
            log.info("[n=%d] PSL size: %.1f KB", n_genomes, psl_file.stat().st_size / 1024)

            # 3. Run GEESE
            geese_file = iter_dir / "atoms.geese"
            run_geese(
                geese_bin,
                psl_file,
                geese_file,
                min_length=min_length,
                min_alignment_length=min_alignment_length,
                min_identity=min_identity,
            )
            log.info("[n=%d] GEESE size: %.1f KB", n_genomes, geese_file.stat().st_size / 1024)

            # 4. Run scorer
            scorer_dir = iter_dir / "scorer_output"
            scores = run_scorer(subset_fasta, geese_file, scorer_dir)

            if scores["overall"] is not None:
                row["overall_score"] = f"{scores['overall']:.6f}"
                row["alignment_score"] = f"{scores['alignment']:.6f}"
                row["coverage_score"] = f"{scores['coverage']:.6f}"
                log.info(
                    "[n=%d] SCORES  overall=%.4f  alignment=%.4f  coverage=%.4f",
                    n_genomes,
                    scores["overall"],
                    scores["alignment"],
                    scores["coverage"],
                )
            else:
                row["error"] = "scorer returned None for one or more sub-scores"
                log.warning("[n=%d] Scorer returned partial/no results", n_genomes)

        except subprocess.CalledProcessError as error:
            stderr_snippet = (error.stderr or "")[:500]
            row["error"] = f"cmd failed (rc={error.returncode}): {stderr_snippet}"
            log.error("[n=%d] FAILED: %s", n_genomes, row["error"])

        except Exception as error:
            row["error"] = str(error)[:500]
            log.error("[n=%d] FAILED: %s", n_genomes, error, exc_info=True)

        results.append(row)
        write_results_csv(results, results_csv)

    log.info("=" * 60)
    log.info("Done.  CSV: %s", results_csv)
    return results_csv


# ---------------------------------------------------------------------------
# Script Entry Point
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-8s  %(message)s",
        datefmt="%H:%M:%S",
    )
    args = parse_args()
    try:
        run(
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
    except ValueError as exc:
        log.error("%s", exc)
        sys.exit(1)
