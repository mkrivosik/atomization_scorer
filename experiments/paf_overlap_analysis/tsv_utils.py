"""Shared TSV writing utility for paf_overlap_analysis modules."""

# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------
from __future__ import annotations

from pathlib import Path


# ---------------------------------------------------------------------------
# TSV Writing
# ---------------------------------------------------------------------------
def write_tsv(output_path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    """Write a padded tab-delimited table with a dashed separator under the header."""
    column_widths = {
        fieldname: max(
            len(fieldname),
            max((len(str(row.get(fieldname, ""))) for row in rows), default=0),
        )
        for fieldname in fieldnames
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as file:
        header = "\t".join(fieldname.ljust(column_widths[fieldname]) for fieldname in fieldnames)
        separator = "\t".join("-" * column_widths[fieldname] for fieldname in fieldnames)
        file.write(header + "\n")
        file.write(separator + "\n")
        for row in rows:
            line = "\t".join(
                str(row.get(fieldname, "")).ljust(column_widths[fieldname])
                for fieldname in fieldnames
            )
            file.write(line + "\n")
