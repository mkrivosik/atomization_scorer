"""Analyze large primary-primary overlaps (1000-9999 bp) from a filtered minimap2 PAF file."""

# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------
from __future__ import annotations

import argparse
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any

from tsv_utils import write_tsv

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
LARGE_PP_MIN = 1000
LARGE_PP_MAX = 9999
RECORDS_FIELDS = [
    "genome_id",
    "class_a",
    "class_b",
    "alignment_start_a",
    "alignment_end_a",
    "alignment_start_b",
    "alignment_end_b",
    "representant_a_length",
    "representant_b_length",
    "overlap_length",
    "overlap_type",
    "strand_pattern",
]
SUMMARY_FIELDS = [
    "class_a",
    "class_b",
    "overlap_type",
    "strand_pattern",
    "length_bin",
    "n_events",
    "n_genomes",
    "median_overlap",
    "sd_overlap",
    "min_overlap",
    "max_overlap",
]


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------
def _parse_class(target_name: str) -> str:
    """Extract class ID from a representative target name like '<genome>|class_<id>'."""
    if "|class_" in target_name:
        return target_name.split("|class_")[1]
    return target_name


def parse_alignments_rich(paf_path: Path) -> list[dict[str, Any]]:
    """Read PAF and return alignment records with query coordinates, strand, target, and type."""
    alignments = []
    with paf_path.open() as file:
        for line in file:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            atype = "unknown"
            for field in fields[12:]:
                if field.startswith("tp:A:"):
                    atype = field.split(":")[2]
                    break
            alignments.append({
                "query_name": fields[0],
                "query_start": int(fields[2]),
                "query_end": int(fields[3]),
                "strand": fields[4],
                "target_name": fields[5],
                "target_length": int(fields[6]),
                "atype": atype,
                "class_id": _parse_class(fields[5]),
            })
    return alignments


# ---------------------------------------------------------------------------
# Overlap Detection
# ---------------------------------------------------------------------------
def _classify_overlap_type(alignment_end_a: int, alignment_end_b: int) -> str:
    """Classify overlap as both_edge (cross-overlap) or mixed_edge_internal (containment)."""
    if alignment_end_a <= alignment_end_b:
        return "both_edge"
    return "mixed_edge_internal"


def detect_large_pp_records(
    alignments: list[dict[str, Any]],
    min_overlap: int = LARGE_PP_MIN,
    max_overlap: int = LARGE_PP_MAX,
) -> list[dict[str, Any]]:
    """Detect P-P pairs in [min_overlap, max_overlap] and return one rich record per pair."""
    grouped = defaultdict(list)
    for alignment in alignments:
        grouped[alignment["query_name"]].append(alignment)

    records = []
    for genome_alignments in grouped.values():
        sorted_alignments = sorted(
            genome_alignments,
            key=lambda aln: (aln["query_start"], aln["query_end"]),
        )
        alignments_count = len(sorted_alignments)
        for first_index in range(alignments_count):
            for second_index in range(first_index + 1, alignments_count):
                alignment = sorted_alignments[first_index]
                partner = sorted_alignments[second_index]
                alignment_start_a = int(alignment["query_start"])
                alignment_end_a = int(alignment["query_end"])
                alignment_start_b = int(partner["query_start"])
                alignment_end_b = int(partner["query_end"])

                if alignment_start_b >= alignment_end_a:
                    break
                if alignment["atype"] != "P" or partner["atype"] != "P":
                    continue

                overlap_length = min(alignment_end_a, alignment_end_b) - max(alignment_start_a, alignment_start_b)
                if overlap_length <= 0:
                    continue
                if not (min_overlap <= overlap_length <= max_overlap):
                    continue

                strand_pattern = "same" if alignment["strand"] == partner["strand"] else "opposite"
                overlap_type = _classify_overlap_type(alignment_end_a, alignment_end_b)

                records.append({
                    "genome_id": alignment["query_name"],
                    "class_a": alignment["class_id"],
                    "class_b": partner["class_id"],
                    "alignment_start_a": alignment_start_a,
                    "alignment_end_a": alignment_end_a,
                    "alignment_start_b": alignment_start_b,
                    "alignment_end_b": alignment_end_b,
                    "representant_a_length": alignment["target_length"],
                    "representant_b_length": partner["target_length"],
                    "overlap_length": overlap_length,
                    "overlap_type": overlap_type,
                    "strand_pattern": strand_pattern,
                })

    return records


# ---------------------------------------------------------------------------
# Signature Building
# ---------------------------------------------------------------------------
def _length_bin(overlap_length: int) -> str:
    """Return a length bin label within the 1000-9999 bp range."""
    if overlap_length < 2000:
        return "1000-1999"
    if overlap_length < 3000:
        return "2000-2999"
    if overlap_length < 4000:
        return "3000-3999"
    if overlap_length < 5000:
        return "4000-4999"
    return "5000-9999"


def build_signatures(records: list[dict[str, Any]]) -> dict[tuple, list[dict[str, Any]]]:
    """Group records by (class_pair, overlap_type, strand_pattern, length_bin)."""
    signatures = defaultdict(list)
    for record in records:
        key = (
            tuple(sorted([record["class_a"], record["class_b"]])),
            record["overlap_type"],
            record["strand_pattern"],
            _length_bin(record["overlap_length"]),
        )
        signatures[key].append(record)
    return dict(signatures)


def summarize_signatures(signatures: dict[tuple, list[dict[str, Any]]]) -> list[dict[str, Any]]:
    """Compute per-signature statistics, sorted by n_genomes then median_overlap descending."""
    rows = []
    for signature_key, records in signatures.items():
        class_pair, overlap_type, strand_pattern, length_bin = signature_key
        class_a, class_b = class_pair
        overlap_lengths = [record["overlap_length"] for record in records]
        rows.append({
            "class_a": class_a,
            "class_b": class_b,
            "overlap_type": overlap_type,
            "strand_pattern": strand_pattern,
            "length_bin": length_bin,
            "n_events": len(records),
            "n_genomes": len({r["genome_id"] for r in records}),
            "median_overlap": round(statistics.median(overlap_lengths)),
            "sd_overlap": round(statistics.stdev(overlap_lengths), 1) if len(records) > 1 else 0.0,
            "min_overlap": min(overlap_lengths),
            "max_overlap": max(overlap_lengths),
        })
    return sorted(rows, key=lambda row: (-row["n_genomes"], -row["median_overlap"]))


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------
def write_outputs(records: list[dict[str, Any]], summaries: list[dict[str, Any]], output_dir: Path) -> None:
    """Write large_pp_records.tsv and large_pp_signatures.tsv to output_dir."""
    records_path = output_dir / "large_pp_records.tsv"
    write_tsv(records_path, records, RECORDS_FIELDS)
    summary_path = output_dir / "large_pp_signatures.tsv"
    write_tsv(summary_path, summaries, SUMMARY_FIELDS)
    print(f"Records    ({len(records)} rows): {records_path}")
    print(f"Signatures ({len(summaries)} rows): {summary_path}")


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
def run(paf_path: Path, output_dir: Path) -> None:
    """Parse PAF, detect large P-P overlaps, group into signatures, write reports."""
    if not paf_path.is_file():
        raise FileNotFoundError(f"PAF file not found: {paf_path}")

    alignments = parse_alignments_rich(paf_path)
    records = detect_large_pp_records(alignments)
    signatures = build_signatures(records)
    summaries = summarize_signatures(signatures)
    write_outputs(records, summaries, output_dir)

    print(f"\nTotal large P-P records: {len(records)}")
    print(f"Distinct signatures:     {len(signatures)}")
    if summaries:
        print("Top signatures by n_genomes:")
        for summary in summaries[:5]:
            print(
                f"  ({summary['class_a']}, {summary['class_b']}) "
                f"{summary['overlap_type']} {summary['strand_pattern']} "
                f"{summary['length_bin']}bp "
                f"-> {summary['n_events']} events in {summary['n_genomes']} genomes"
            )


# ---------------------------------------------------------------------------
# Command Line Interface
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Analyze large P-P overlaps (1000-9999 bp) and group into signatures."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("input/minimap2_alignment_filtered.paf"),
        help="Path to the filtered PAF file (default: input/minimap2_alignment_filtered.paf)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("reports/large_pp"),
        help="Output directory for TSV reports (default: reports/large_pp)",
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Script Entry Point
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    run(args.input, args.output)
