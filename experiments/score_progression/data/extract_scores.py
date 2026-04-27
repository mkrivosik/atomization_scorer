"""
Extract the four score columns from the results CSV and write an aligned scores.tsv.
"""

# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------
from __future__ import annotations

import argparse
import csv
from pathlib import Path

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
COLUMNS = ["n_genomes", "overall_score", "alignment_score", "coverage_score"]


# ---------------------------------------------------------------------------
# Command Line Interface
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Extract score columns from the results CSV into an aligned TSV.")
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Path to the source results CSV",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Path for the output TSV (default: scores.tsv next to --input)",
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
def run(input_csv: Path, output_tsv: Path) -> Path:
    """Extract score columns from the results CSV and write an aligned scores.tsv."""
    if not input_csv.is_file():
        raise FileNotFoundError(f"Input CSV not found: {input_csv}")

    with open(input_csv, newline="") as file:
        rows = list(csv.DictReader(file))

    kept: list[dict[str, str]] = []
    for row in rows:
        if row.get("error"):
            continue
        if not all(row.get(column) for column in COLUMNS):
            continue
        kept.append({column: row[column] for column in COLUMNS})

    if not kept:
        print("No valid rows found — scores.tsv not written.")
        return output_tsv

    widths = {column: len(column) for column in COLUMNS}
    for row in kept:
        for column in COLUMNS:
            widths[column] = max(widths[column], len(row[column]))

    def pad_row(values: list[str]) -> str:
        return "\t".join(value.ljust(widths[column]) for column, value in zip(COLUMNS, values)).rstrip()

    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(output_tsv, "w") as file:
        file.write(pad_row(COLUMNS) + "\n")
        for row in kept:
            file.write(pad_row([row[column] for column in COLUMNS]) + "\n")

    print(f"Written {len(kept)} rows to {output_tsv}")
    return output_tsv


# ---------------------------------------------------------------------------
# Script Entry Point
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    output = args.output or args.input.parent / "scores.tsv"
    run(args.input, output)
