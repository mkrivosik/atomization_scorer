"""
Create a subset of a FASTA file by removing sequences at specified 1-based positions.
"""

# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------
from __future__ import annotations

import argparse
import logging
from pathlib import Path

from Bio import SeqIO

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
_HERE                   = Path(__file__).resolve().parent
FASTA                   = _HERE / "../../../tests/fixtures/big.fa"
OUTPUT_DIR              = _HERE / "../../../tests/fixtures"
# Score progression on big.fa (598 genomes) revealed three additions that caused sharp drops:
#   position 230 -> SAMN15148346-u.1 (chromosome): score 0.933 -> 0.768
#   position 232 -> SAMN15148347-u.1 (chromosome): score 0.768 -> 0.608
#   position 590 -> SAMN15148663-u.1 (chromosome): score 0.584 -> 0.456
#
# A fourth drop was found by running score progression on good_genomes.fa (big.fa minus the
# three above). Because removing positions 230, 232, 590 shifts every subsequent index by -3,
# a drop observed at position 591 in good_genomes.fa maps back to position 591 + 3 = 594 in big.fa.
#   position 594 -> SAMN15148693-u.1 (chromosome): score 0.922 -> 0.563
# Removing all four outlier chromosomes and scoring good_genomes.fa produces a stable score of ~0.92,
# indicating high-quality atomization on the remaining genomes. This was further validated by an
# artificial degradation test: randomly reassigning the class of x*10% of atoms (x = 1...10)
# in the predicted GEESE consistently lowered the score, with larger fractions producing larger
# drops, confirming the metric penalises class misassignments and that the ~0.92 baseline
# reflects real atomization quality.
DEFAULT_EXCLUDE_INDICES = [230, 232, 590, 594]


# ---------------------------------------------------------------------------
# Create Subset FASTA
# ---------------------------------------------------------------------------
def create_subset_fasta(
    input_fasta: Path,
    output_dir: Path,
    exclude_indices: list[int] | None = None,
) -> Path:
    """Remove sequences by 1-based position and write the result to good_genomes.fa."""
    output_dir.mkdir(parents=True, exist_ok=True)
    resolved = exclude_indices if exclude_indices is not None else DEFAULT_EXCLUDE_INDICES

    records = list(SeqIO.parse(input_fasta, "fasta"))
    exclude_set = set(resolved)

    kept = []
    deleted = []
    for index, record in enumerate(records, start=1):
        if index in exclude_set:
            deleted.append(record)
        else:
            kept.append(record)

    output_path = output_dir / "good_genomes.fa"
    SeqIO.write(kept, output_path, "fasta")
    log.info("Wrote %d sequences (excluded %d) -> %s", len(kept), len(deleted), output_path)

    if deleted:
        deleted_path = output_dir / "bad_genomes.fa"
        SeqIO.write(deleted, deleted_path, "fasta")
        log.info("Wrote %d deleted sequences -> %s", len(deleted), deleted_path)

    return output_path


# ---------------------------------------------------------------------------
# Script Entry Point
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-8s  %(message)s",
        datefmt="%H:%M:%S",
    )
    parser = argparse.ArgumentParser(
        description="Create a subset FASTA by dropping sequences at given 1-based positions."
    )
    parser.add_argument(
        "input_fasta",
        type=Path,
        nargs="?",
        default=FASTA,
        help=f"Source FASTA file (default: {FASTA})",
    )
    parser.add_argument(
        "output_dir",
        type=Path,
        nargs="?",
        default=OUTPUT_DIR,
        help=f"Output directory (default: {OUTPUT_DIR})",
    )
    parser.add_argument(
        "--exclude",
        type=int,
        nargs="+",
        default=DEFAULT_EXCLUDE_INDICES,
        metavar="N",
        help=f"1-based sequence indices to exclude (default: {DEFAULT_EXCLUDE_INDICES})",
    )
    args = parser.parse_args()

    create_subset_fasta(args.input_fasta, args.output_dir, args.exclude)
