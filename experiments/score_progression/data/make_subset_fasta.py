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
_HERE                   = Path(__file__).parent
FASTA                   = _HERE / "../../tests/fixtures/big.fa"
OUTPUT_DIR              = _HERE / "../../tests/fixtures"
DEFAULT_EXCLUDE_INDICES = [] # 230, 232


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
