"""
Parse a filtered minimap2 PAF file and write a JSON overlap report broken down by
pair type (P-P, P-S, S-S) and overlap length category.
"""

# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
PAIR_TYPES = ["P-P", "P-S", "S-S"]
LENGTH_CATEGORIES = [
    ("1_to_9bp",        1,      9),
    ("10_to_99bp",      10,     99),
    ("100_to_999bp",    100,    999),
    ("1000_to_9999bp",  1000,   9999),
    (">=10000bp",       10000,  None),
]


# ---------------------------------------------------------------------------
# Command Line Interface
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Parse overlapping PAF alignment pairs and write a JSON report."
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
        default=Path("reports/overlap_report.json"),
        help="Path for the output JSON report (default: reports/overlap_report.json)",
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------
def parse_alignments(paf_path: Path) -> list[dict]:
    """Read PAF rows and return a list of alignment records with query coordinates and type."""
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
                "atype": atype,
            })
    return alignments


# ---------------------------------------------------------------------------
# Overlap Detection
# ---------------------------------------------------------------------------
def detect_overlapping_pairs(alignments: list[dict]) -> list[tuple[int, str]]:
    """Return (overlap_length, pair_type) for every overlapping alignment pair."""
    grouped = defaultdict(list)
    for alignment in alignments:
        grouped[alignment["query_name"]].append(alignment)

    pairs = []
    for alignments in grouped.values():
        sorted_alignments = sorted(alignments, key=lambda aln: (aln["query_start"], aln["query_end"]))
        alignments_count = len(sorted_alignments)
        for first_index in range(alignments_count):
            for second_index in range(first_index + 1, alignments_count):
                alignment = sorted_alignments[first_index]
                partner = sorted_alignments[second_index]
                if partner["query_start"] >= alignment["query_end"]:
                    break
                overlap = (
                    min(int(alignment["query_end"]), int(partner["query_end"]))
                    - max(int(alignment["query_start"]), int(partner["query_start"]))
                )
                if overlap <= 0:
                    continue
                first_type = str(alignment["atype"])
                second_type = str(partner["atype"])
                if first_type == "P" and second_type == "P":
                    pair_type = "P-P"
                elif first_type == "S" and second_type == "S":
                    pair_type = "S-S"
                elif (first_type == "P" and second_type == "S") or (first_type == "S" and second_type == "P"):
                    pair_type = "P-S"
                else:
                    pair_type = "other"
                pairs.append((overlap, pair_type))

    return pairs


# ---------------------------------------------------------------------------
# Report Building
# ---------------------------------------------------------------------------
def _length_category(length: int) -> str:
    """Return the length-category label for a given overlap length."""
    for name, low, high in LENGTH_CATEGORIES:
        if high is None or length <= high:
            return name
    return ">=10000bp"


def build_report(pairs: list[tuple[int, str]], total_alignments: int) -> dict:
    """Aggregate overlapping pairs into a nested report dictionary."""
    empty_categories = {name: 0 for name, *_ in LENGTH_CATEGORIES}
    counts = {pair_type: {"total": 0, **empty_categories.copy()} for pair_type in PAIR_TYPES}
    for overlap_length, pair_type in pairs:
        if pair_type not in counts:
            continue
        counts[pair_type]["total"] += 1
        counts[pair_type][_length_category(overlap_length)] += 1
    return {
        "total_alignments": total_alignments,
        "total_overlapping_pairs": len(pairs),
        "P-P": counts["P-P"],
        "P-S": counts["P-S"],
        "S-S": counts["S-S"],
    }


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
def run(paf_path: Path, output_path: Path) -> None:
    """Parse PAF overlaps and write the JSON report to output_path."""
    if not paf_path.is_file():
        raise FileNotFoundError(f"PAF file not found: {paf_path}")

    alignments = parse_alignments(paf_path)
    pairs = detect_overlapping_pairs(alignments)
    report = build_report(pairs, total_alignments=len(alignments))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as file:
        json.dump(report, file, indent=2)
        file.write("\n")

    print(f"Total alignments:        {report['total_alignments']}")
    print(f"Total overlapping pairs: {report['total_overlapping_pairs']}")
    for pair_type in PAIR_TYPES:
        print(f"  {pair_type}: {report[pair_type]['total']}")
    print(f"Report written to {output_path}")


# ---------------------------------------------------------------------------
# Script Entry Point
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    run(args.input, args.output)
