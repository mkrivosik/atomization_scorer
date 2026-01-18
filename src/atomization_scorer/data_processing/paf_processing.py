"""
paf_processing.py

Utility functions for processing Minimap2 PAF outputs.

Functions
---------
filter_paf : Filter PAF alignments by minimum similarity and minimum alignment length.
"""

# ---------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------
from pathlib import Path

# ---------------------------------------------------------------------
# PAF Processing
# ---------------------------------------------------------------------

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
        Path to input PAF file from minimap2.
    output_file : Path
        Path to the output directory where filtered PAF file is stored.
    minimum_similarity : float, optional, default=0.95
        Minimum fraction of matching bases (0.0 to 1.0).
    minimum_alignment_length : int, optional, default=500
        Minimum number of aligned bases.

    Raises
    ------
    FileNotFoundError
        Raised if the input PAF file does not exist.

    Returns
    -------
    Path
        Path to filtered PAF file.
    """
    if not paf_file.is_file():
        raise FileNotFoundError(f"PAF file not found: {paf_file}")

    output_file.parent.mkdir(parents=True, exist_ok=True)

    filtered_lines = []

    with paf_file.open("r") as file:
        for line in file:
            fields = line.strip().split("\t")

            alignment_length = int(fields[10])

            # Divergence: find optional field starting with "de:f:"
            for field in fields[12:]:
                if field.startswith("de:f:"):
                    divergence = float(field.split(":")[2])
                    break

            similarity = 1.0 - divergence

            if similarity >= minimum_similarity and alignment_length >= minimum_alignment_length:
                filtered_lines.append(line)

    with output_file.open("w") as out:
        out.writelines(filtered_lines)

    return output_file
