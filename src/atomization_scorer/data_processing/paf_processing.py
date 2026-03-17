"""
paf_processing.py

Utility functions for processing Minimap2 PAF outputs.

Functions
---------
filter_paf : Filter PAF alignments by minimum similarity and minimum alignment length.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
import logging
from pathlib import Path

# --------------------------------------------------------------------------------------
# Logging
# --------------------------------------------------------------------------------------
logger = logging.getLogger(__name__)

# --------------------------------------------------------------------------------------
# PAF Processing
# --------------------------------------------------------------------------------------

def filter_paf(
    paf_file: Path,
    output_file: Path,
    minimum_similarity: float = 0.95,
    minimum_alignment_length: int = 500
) -> Path:
    """
    Filter PAF alignments based on minimum similarity and minimum alignment length.

    Parameters
    ----------
    paf_file : Path
        Path to the input PAF file from minimap2.
    output_file : Path
        Path to the filtered PAF file.
    minimum_similarity : float, optional, default=0.95
        Minimum fraction of matching bases (0.0 to 1.0).
    minimum_alignment_length : int, optional, default=500
        Minimum number of aligned bases.

    Raises
    ------
    FileNotFoundError
        Raised if the input PAF file does not exist.
    ValueError
        Raised if the PAF file contains malformed lines or invalid numeric fields.

    Returns
    -------
    Path
        Filtered PAF file path.
    """
    if not paf_file.is_file():
        raise FileNotFoundError(f"PAF file not found: {paf_file}")

    output_file.parent.mkdir(parents=True, exist_ok=True)

    filtered_lines = []
    total_lines = 0

    logger.info(
        "Filtering PAF file %s into %s with minimum_similarity=%s and minimum_alignment_length=%s",
        paf_file,
        output_file,
        minimum_similarity,
        minimum_alignment_length,
    )

    with paf_file.open("r") as file:
        for line in file:
            total_lines += 1
            fields = line.strip().split("\t")
            if len(fields) < 11:
                raise ValueError(f"Malformed PAF line: {line.strip()}")

            try:
                alignment_length = int(fields[10])
                matches = int(fields[9])
            except ValueError as error:
                raise ValueError(f"Malformed PAF line: {line.strip()}") from error

            # Prefer minimap2's dv tag; fall back to the identity estimate from mandatory PAF columns.
            similarity = matches / alignment_length if alignment_length > 0 else 0.0

            for field in fields[12:]:
                if field.startswith("dv:f:"):
                    similarity = 1.0 - float(field.split(":")[2])
                    break

            if similarity >= minimum_similarity and alignment_length >= minimum_alignment_length:
                filtered_lines.append(line)

    with output_file.open("w") as out:
        out.writelines(filtered_lines)

    logger.info(
        "Filtered PAF saved to %s with %s kept alignments out of %s total alignments",
        output_file,
        len(filtered_lines),
        total_lines,
    )

    return output_file
