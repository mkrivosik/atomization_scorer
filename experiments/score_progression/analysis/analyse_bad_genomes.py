"""
Analyse atom statistics for the two bad genomes identified in score_progression experiments.
"""

# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------
from __future__ import annotations

import argparse
import json
import logging
import statistics
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from atomization_scorer.data_processing import read_geese

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
HERE = Path(__file__).parent
FIXTURES = HERE.parent.parent.parent / "tests" / "fixtures"

DEFAULT_FASTA = FIXTURES / "mini.fa"
DEFAULT_GEESE = FIXTURES / "mini.geese"
DEFAULT_DELETED = FIXTURES / "deleted_v1.fa"
DEFAULT_OUTPUT = HERE.parent / "outputs" / "bad_genomes_report.json"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _median(values: list[float]) -> float | None:
    """Return median of values, or None if empty."""
    return statistics.median(values) if values else None


def _mean(values: list[float]) -> float | None:
    """Return mean of values, or None if empty."""
    return statistics.mean(values) if values else None


def _percentile(values: list[float], p: float) -> float | None:
    """Return the p-th percentile of values using linear interpolation, or None if empty."""
    if not values:
        return None
    sorted_values = sorted(values)
    index = (p / 100) * (len(sorted_values) - 1)
    lower = int(index)
    upper = lower + 1
    if upper >= len(sorted_values):
        return float(sorted_values[-1])
    fraction = index - lower
    return float(sorted_values[lower] + fraction * (sorted_values[upper] - sorted_values[lower]))


def _compute_reference_atom_length(atoms_df: pd.DataFrame) -> float:
    """Return reference atom length: median of good atoms (multi-class globally, length ≤ global P99)."""
    df = atoms_df.copy()
    df["length"] = df["end"].astype(int) - df["start"].astype(int)

    class_counts = df.groupby("class").size()
    multi_classes = set(class_counts[class_counts >= 2].index)
    good_lengths = df[df["class"].isin(multi_classes)]["length"].tolist()

    p99 = _percentile(good_lengths, 99)
    if p99 is None:
        return statistics.median(good_lengths)
    filtered = [length for length in good_lengths if length <= p99]
    return statistics.median(filtered)


def _merge_intervals(intervals: list[tuple[int, int]]) -> int:
    """Return total covered length after merging overlapping half-open intervals."""
    if not intervals:
        return 0
    sorted_intervals = sorted(intervals)
    current_start, current_end = sorted_intervals[0]
    total = 0
    for start, end in sorted_intervals[1:]:
        if start <= current_end:
            current_end = max(current_end, end)
        else:
            total += current_end - current_start
            current_start, current_end = start, end
    total += current_end - current_start
    return total


# ---------------------------------------------------------------------------
# PAF Parsing
# ---------------------------------------------------------------------------
def _parse_paf(paf_file: Path) -> dict[str, list[tuple[int, int, str]]]:
    """Parse PAF into per-genome alignment records: {genome: [(start, end, class_id), ...]}."""
    alignments = defaultdict(list)
    with paf_file.open() as file:
        for line in file:
            if not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) < 12:
                continue
            query_name = fields[0]
            query_start = int(fields[2])
            query_end = int(fields[3])
            target_name = fields[5]
            if "|class_" not in target_name:
                continue
            class_id = target_name.split("|class_")[1].strip()
            alignments[query_name].append((query_start, query_end, class_id))
    return dict(alignments)


# ---------------------------------------------------------------------------
# Atom Statistics
# ---------------------------------------------------------------------------
def _compute_large_atom_pcts(lengths: pd.Series, total_atoms: int) -> tuple[float, float, float]:
    """Return (pct_above_10kb, pct_above_50kb, pct_above_100kb)."""
    if total_atoms == 0:
        return 0.0, 0.0, 0.0
    return (
        round(int((lengths > 10_000).sum()) / total_atoms * 100, 4),
        round(int((lengths > 50_000).sum()) / total_atoms * 100, 4),
        round(int((lengths > 100_000).sum()) / total_atoms * 100, 4),
    )


def _compute_alignment_support_stats(
    atoms_df: pd.DataFrame,
    paf_alignments: dict[str, list[tuple[int, int, str]]],
    global_singleton_str: set[str],
    total_atoms: int,
) -> dict:
    """Count alignment support metrics for a set of atoms and return stats dictionary."""
    n_no_support = 0
    n_strong_support = 0
    n_singleton_total = 0
    n_singleton_with_support = 0

    for _, row in atoms_df.iterrows():
        genome_alignments = paf_alignments.get(row["name"], [])
        n_alignments, _ = _compute_atom_alignment_support(
            atom_start=int(row["start"]),
            atom_end=int(row["end"]),
            atom_class=str(row["class"]),
            genome_alignments=genome_alignments,
        )
        if n_alignments == 0:
            n_no_support += 1
        if n_alignments >= 3:
            n_strong_support += 1
        if str(row["class"]) in global_singleton_str:
            n_singleton_total += 1
            if n_alignments >= 1:
                n_singleton_with_support += 1

    return {
        "pct_atoms_no_alignment_support": round(n_no_support / total_atoms * 100, 4) if total_atoms > 0 else 0.0,
        "pct_atoms_strong_support": round(n_strong_support / total_atoms * 100, 4) if total_atoms > 0 else 0.0,
        "pct_singleton_atoms_with_alignment_support": round(n_singleton_with_support / n_singleton_total * 100, 4) if n_singleton_total > 0 else 0.0,
    }


def _compute_atom_alignment_support(
    atom_start: int,
    atom_end: int,
    atom_class: str,
    genome_alignments: list[tuple[int, int, str]],
) -> tuple[int, float]:
    """Return (n_alignments, coverage_fraction) for one atom."""
    atom_length = atom_end - atom_start
    if atom_length == 0:
        return 0, 0.0

    overlapping = [
        (max(atom_start, alignment_start), min(atom_end, alignment_end))
        for alignment_start, alignment_end, alignment_class in genome_alignments
        if alignment_class == atom_class and min(atom_end, alignment_end) > max(atom_start, alignment_start)
    ]

    n_alignments = len(overlapping)
    covered = _merge_intervals(overlapping)
    return n_alignments, covered / atom_length


# ---------------------------------------------------------------------------
# Per-Genome Analysis
# ---------------------------------------------------------------------------
def _analyse_genome(
    genome_name: str,
    genome_length: int,
    atoms_df: pd.DataFrame,
    global_singleton_classes: set,
    reference_atom_length: float,
    paf_alignments: dict[str, list[tuple[int, int, str]]] | None,
) -> dict:
    """Compute all stats for one bad genome."""
    genome_atoms = atoms_df[atoms_df["name"] == genome_name].copy()
    genome_atoms["length"] = genome_atoms["end"].astype(int) - genome_atoms["start"].astype(int)

    class_counts = genome_atoms.groupby("class").size()
    single_classes = set(class_counts[class_counts == 1].index)
    multi_classes = set(class_counts[class_counts >= 2].index)

    single_atoms = genome_atoms[genome_atoms["class"].isin(single_classes)]
    multi_atoms = genome_atoms[genome_atoms["class"].isin(multi_classes)]

    single_lengths = single_atoms["length"].tolist()
    multi_lengths = multi_atoms["length"].tolist()
    all_lengths = genome_atoms["length"].tolist()

    covered_bp = int(genome_atoms["length"].sum())
    single_bp = int(single_atoms["length"].sum())
    multi_bp = int(multi_atoms["length"].sum())
    uncovered_bp = genome_length - covered_bp
    total_atoms = len(genome_atoms)

    avg_atoms_per_multi = (
        statistics.mean(class_counts[class_counts >= 2].tolist())
        if multi_classes else None
    )

    # Fragmentation score
    expected_atoms = genome_length / reference_atom_length
    fragmentation_score = round(total_atoms / expected_atoms, 4) if expected_atoms > 0 else None

    # Global singleton proportion
    n_in_global_singletons = int(genome_atoms["class"].isin(global_singleton_classes).sum())
    pct_in_global_singletons = round(n_in_global_singletons / total_atoms * 100, 4) if total_atoms > 0 else 0.0

    # Large atom proportions
    pct_above_10kb, pct_above_50kb, pct_above_100kb = _compute_large_atom_pcts(genome_atoms["length"], total_atoms)

    # Length percentiles
    median_len = _median(all_lengths)
    p90 = _percentile(all_lengths, 90)
    p95 = _percentile(all_lengths, 95)
    p99 = _percentile(all_lengths, 99)
    large_atom_warning = bool(p99 is not None and median_len and p99 / median_len >= 10)

    result = {
        "genome_name": genome_name,
        "genome_length_bp": genome_length,

        "total_atoms": total_atoms,
        "atoms_in_single_atom_classes": len(single_atoms),
        "atoms_in_multi_atom_classes": len(multi_atoms),

        "total_distinct_classes": len(class_counts),
        "n_single_atom_classes": len(single_classes),
        "n_multi_atom_classes": len(multi_classes),

        "covered_bp": covered_bp,
        "single_atom_class_bp": single_bp,
        "multi_atom_class_bp": multi_bp,
        "uncovered_bp": uncovered_bp,

        "coverage_pct": round(covered_bp / genome_length * 100, 4),
        "single_atom_class_pct": round(single_bp / genome_length * 100, 4),
        "multi_atom_class_pct": round(multi_bp / genome_length * 100, 4),
        "uncovered_pct": round(uncovered_bp / genome_length * 100, 4),

        "median_atom_length": median_len,
        "p90_atom_length": p90,
        "p95_atom_length": p95,
        "p99_atom_length": p99,
        "large_atom_warning": large_atom_warning,
        "median_single_atom_class_length": _median(single_lengths),
        "median_multi_atom_class_atom_length": _median(multi_lengths),
        "avg_atoms_per_multi_class": round(avg_atoms_per_multi, 4) if avg_atoms_per_multi is not None else None,

        "expected_atoms": round(expected_atoms, 2),
        "fragmentation_score": fragmentation_score,

        "pct_atoms_in_global_singleton_classes": pct_in_global_singletons,

        "pct_atoms_above_10kb": pct_above_10kb,
        "pct_atoms_above_50kb": pct_above_50kb,
        "pct_atoms_above_100kb": pct_above_100kb,
    }

    if paf_alignments is not None:
        global_singleton_str = {str(c) for c in global_singleton_classes}
        result.update(_compute_alignment_support_stats(
            atoms_df=genome_atoms,
            paf_alignments=paf_alignments,
            global_singleton_str=global_singleton_str,
            total_atoms=total_atoms,
        ))

    return result


# ---------------------------------------------------------------------------
# Summary across all genomes
# ---------------------------------------------------------------------------
def _build_summary(
    atoms_df: pd.DataFrame,
    paf_alignments: dict[str, list[tuple[int, int, str]]] | None,
    global_singleton_classes: set,
    reference_atom_length: float,
) -> dict:
    """Compute aggregate length stats across all genomes in the GEESE file."""
    all_atoms = atoms_df.copy()
    all_atoms["length"] = all_atoms["end"].astype(int) - all_atoms["start"].astype(int)

    all_lengths = all_atoms["length"].tolist()
    single_lengths = []
    multi_lengths = []
    multi_class_sizes = []

    for genome_name in all_atoms["name"].unique():
        genome_atoms = all_atoms[all_atoms["name"] == genome_name]
        class_counts = genome_atoms.groupby("class").size()
        single_classes = set(class_counts[class_counts == 1].index)
        multi_classes = set(class_counts[class_counts >= 2].index)
        single_lengths.extend(genome_atoms[genome_atoms["class"].isin(single_classes)]["length"].tolist())
        multi_lengths.extend(genome_atoms[genome_atoms["class"].isin(multi_classes)]["length"].tolist())
        multi_class_sizes.extend(class_counts[class_counts >= 2].tolist())

    # Class sharing
    class_sharing = atoms_df.groupby("class")["name"].nunique().to_dict()
    genomes_per_class = list(class_sharing.values())
    n_total_classes = len(genomes_per_class)
    n_singleton_classes = int((atoms_df.groupby("class").size() == 1).sum())
    n_genome_unique_classes = sum(1 for value in genomes_per_class if value == 1)

    # Large atom proportions
    total_atoms = len(all_atoms)
    pct_above_10kb, pct_above_50kb, pct_above_100kb = _compute_large_atom_pcts(all_atoms["length"], total_atoms)

    # Length percentiles
    median_len = _median(all_lengths)
    p99 = _percentile(all_lengths, 99)
    large_atom_warning = bool(p99 is not None and median_len and p99 / median_len >= 10)

    summary = {
        "avg_atom_length": round(_mean(all_lengths) or 0.0, 4) if all_lengths else None,
        "median_atom_length": median_len,
        "p90_atom_length": _percentile(all_lengths, 90),
        "p95_atom_length": _percentile(all_lengths, 95),
        "p99_atom_length": p99,
        "large_atom_warning": large_atom_warning,
        "avg_single_atom_class_length": round(_mean(single_lengths) or 0.0, 4) if single_lengths else None,
        "median_single_atom_class_length": _median(single_lengths),
        "avg_multi_atom_class_atom_length": round(_mean(multi_lengths) or 0.0, 4) if multi_lengths else None,
        "median_multi_atom_class_atom_length": _median(multi_lengths),
        "avg_atoms_per_multi_class": round(statistics.mean(multi_class_sizes), 4) if multi_class_sizes else None,
        "min_atom_length": int(min(all_lengths)) if all_lengths else None,
        "max_atom_length": int(max(all_lengths)) if all_lengths else None,

        "avg_genomes_per_class": round(_mean(genomes_per_class) or 0.0, 4) if genomes_per_class else None,
        "median_genomes_per_class": _median(genomes_per_class),
        "p90_genomes_per_class": _percentile(genomes_per_class, 90),
        "n_total_classes": n_total_classes,
        "n_singleton_classes": n_singleton_classes,
        "pct_singleton_classes": round(n_singleton_classes / n_total_classes * 100, 4) if n_total_classes > 0 else 0.0,
        "n_genome_unique_classes": n_genome_unique_classes,
        "pct_genome_unique_classes": round(n_genome_unique_classes / n_total_classes * 100, 4) if n_total_classes > 0 else 0.0,

        "pct_atoms_above_10kb": pct_above_10kb,
        "pct_atoms_above_50kb": pct_above_50kb,
        "pct_atoms_above_100kb": pct_above_100kb,

        "reference_atom_length_bp": round(reference_atom_length, 4),
    }

    if paf_alignments is not None:
        global_singleton_str = {str(c) for c in global_singleton_classes}
        summary.update(_compute_alignment_support_stats(
            atoms_df=all_atoms,
            paf_alignments=paf_alignments,
            global_singleton_str=global_singleton_str,
            total_atoms=total_atoms,
        ))

    return summary


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
def run(
    fasta: Path,
    geese: Path,
    deleted: Path,
    output: Path,
    paf: Path | None = None,
) -> None:
    """Run bad-genome analysis and write JSON report."""
    log.info("Reading genome lengths from %s", fasta)
    genome_lengths = {
        record.id: len(record.seq)
        for record in SeqIO.parse(fasta, "fasta")
    }

    log.info("Reading bad genome names from %s", deleted)
    bad_genome_names = [record.id for record in SeqIO.parse(deleted, "fasta")]
    log.info("Bad genomes: %s", bad_genome_names)

    log.info("Reading atoms from %s", geese)
    atoms_df = read_geese(geese_file=geese)
    atoms_df["start"] = atoms_df["start"].astype(int)
    atoms_df["end"] = atoms_df["end"].astype(int)

    paf_alignments = None
    if paf is not None:
        log.info("Parsing PAF alignments from %s", paf)
        paf_alignments = _parse_paf(paf)

    class_sharing = atoms_df.groupby("class")["name"].nunique().to_dict()
    global_singleton_classes = {cls for cls, n in class_sharing.items() if n == 1}
    log.info("Global singleton classes: %s of %s total", len(global_singleton_classes), len(class_sharing))

    reference_atom_length = _compute_reference_atom_length(atoms_df)
    log.info("Reference atom length: %.1f bp", reference_atom_length)

    genome_results = []
    for genome_name in bad_genome_names:
        if genome_name not in genome_lengths:
            raise ValueError(f"Genome '{genome_name}' not found in FASTA: {fasta}")
        log.info("Analysing genome: %s", genome_name)
        result = _analyse_genome(
            genome_name=genome_name,
            genome_length=genome_lengths[genome_name],
            atoms_df=atoms_df,
            global_singleton_classes=global_singleton_classes,
            reference_atom_length=reference_atom_length,
            paf_alignments=paf_alignments,
        )
        genome_results.append(result)

    summary = _build_summary(
        atoms_df=atoms_df,
        paf_alignments=paf_alignments,
        global_singleton_classes=global_singleton_classes,
        reference_atom_length=reference_atom_length,
    )

    report = {"genomes": genome_results, "summary": summary}

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w") as file:
        json.dump(report, file, indent=2)
        file.write("\n")

    log.info("Report written to %s", output)


# ---------------------------------------------------------------------------
# Command Line Interface
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Analyse bad genomes from score_progression experiments.")
    parser.add_argument("--fasta", type=Path, default=DEFAULT_FASTA, help="Genomes FASTA file.")
    parser.add_argument("--geese", type=Path, default=DEFAULT_GEESE, help="Atoms GEESE file.")
    parser.add_argument("--deleted", type=Path, default=DEFAULT_DELETED, help="Deleted genomes FASTA file.")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="Output JSON report path.")
    parser.add_argument("--paf", type=Path, default=None, help="Filtered PAF file for alignment support analysis.")
    return parser.parse_args()


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
    run(fasta=args.fasta, geese=args.geese, deleted=args.deleted, output=args.output, paf=args.paf)
