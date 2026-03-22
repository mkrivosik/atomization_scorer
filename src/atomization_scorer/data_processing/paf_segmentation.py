"""
paf_segmentation.py

Resolve filtered PAF alignments into non-overlapping query segmentation.

Functions
---------
_parse_paf_alignment            : Parse one PAF line into a validated alignment record.
_rank_alignment                 : Build a sorting key for prioritizing candidate alignments.
_overlaps_existing              : Check whether an alignment overlaps already accepted query segments.
resolve_paf_overlaps            : Resolve overlapping filtered PAF alignments into valid
                                  non-overlapping query segmentation.
validate_non_overlapping_geese  : Validate that a GEESE file contains non-overlapping
                                  half-open intervals per genome.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
from __future__ import annotations
from dataclasses import dataclass
import logging
from pathlib import Path

from .geese_reader import read_geese

# --------------------------------------------------------------------------------------
# Logging
# --------------------------------------------------------------------------------------
logger = logging.getLogger(__name__)


# --------------------------------------------------------------------------------------
# PAF Alignment Record
# --------------------------------------------------------------------------------------
@dataclass(frozen=True)
class PafAlignment:
    """
    Store parsed PAF alignment fields required for query-segmentation ranking.

    Attributes
    ----------
    fields : List[str]
        Original tab-separated PAF fields.
    query_name : str
        Name of the query genome sequence.
    query_start : int
        Start position of the query interval.
    query_end : int
        End position of the query interval.
    target_name : str
        Name of the target representative sequence.
    matches : int
        Number of matching bases in the alignment.
    alignment_length : int
        Length of the aligned block.
    mapq : int
        Minimap2 mapping quality score.
    similarity : float
        Estimated alignment similarity.
    alignment_type : str
        Alignment type from the minimap2 tp tag.
    """
    fields: list[str]
    query_name: str
    query_start: int
    query_end: int
    target_name: str
    matches: int
    alignment_length: int
    mapq: int
    similarity: float
    alignment_type: str


# --------------------------------------------------------------------------------------
# PAF Alignment Parsing
# --------------------------------------------------------------------------------------
def _parse_paf_alignment(line: str, line_number: int) -> PafAlignment:
    """
    Parse a PAF row into a structured alignment record for query-segmentation ranking.

    Parameters
    ----------
    line : str
        One tab-separated PAF row.
    line_number : int
        Line number of the PAF row in the input file.

    Raises
    ------
    ValueError
        Raised if the PAF row is malformed or contains invalid numeric values.

    Returns
    -------
    PafAlignment
        Parsed PAF alignment record with ranking fields.
    """
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 12:
        raise ValueError(
            f"Malformed PAF row at line {line_number}: expected at least 12 tab-separated fields."
        )

    try:
        query_start = int(fields[2])
        query_end = int(fields[3])
        matches = int(fields[9])
        alignment_length = int(fields[10])
        mapq = int(fields[11])
    except ValueError as error:
        raise ValueError(
            f"Malformed PAF row at line {line_number}: invalid numeric value in mandatory fields."
        ) from error

    if query_start >= query_end:
        raise ValueError(
            f"Malformed PAF row at line {line_number}: query start must be smaller than query end."
        )

    similarity = matches / alignment_length if alignment_length > 0 else 0.0
    alignment_type = "unknown"
    for field in fields[12:]:
        if field.startswith("dv:f:"):
            similarity = 1.0 - float(field.split(":")[2])
        elif field.startswith("tp:A:"):
            alignment_type = field.split(":")[2]

    return PafAlignment(
        fields=fields,
        query_name=fields[0],
        query_start=query_start,
        query_end=query_end,
        target_name=fields[5],
        matches=matches,
        alignment_length=alignment_length,
        mapq=mapq,
        similarity=similarity,
        alignment_type=alignment_type,
    )


# --------------------------------------------------------------------------------------
# PAF Alignment Ranking
# --------------------------------------------------------------------------------------
def _rank_alignment(alignment: PafAlignment) -> tuple[int, int, float, int, int, str, int, int]:
    """
    Build a deterministic ranking tuple for comparing candidate PAF alignments.

    Parameters
    ----------
    alignment : PafAlignment
        Parsed PAF alignment record.

    Returns
    -------
    tuple[int, int, float, int, int, str, int, int]
        Ranking tuple ordered by alignment type priority, mapping quality, similarity,
        alignment length, number of matches, and deterministic tie-breakers.
    """
    alignment_priority = {
        "P": 0,  # primary
        "S": 1,  # secondary
        "I": 2,  # inversion
        "unknown": 3,
    }
    return (
        alignment_priority.get(alignment.alignment_type, 4),
        -alignment.mapq,
        -alignment.similarity,
        -alignment.alignment_length,
        -alignment.matches,
        alignment.target_name,
        alignment.query_start,
        alignment.query_end,
    )


# --------------------------------------------------------------------------------------
# Interval Overlap Check
# --------------------------------------------------------------------------------------
def _overlaps_existing(
    query_start: int,
    query_end: int,
    accepted_intervals: list[tuple[int, int]],
) -> bool:
    """
    Check whether a half-open query interval overlaps any previously accepted interval.

    Parameters
    ----------
    query_start : int
        Start position of the candidate query interval.
    query_end : int
        End position of the candidate query interval.
    accepted_intervals : List[Tuple[int, int]]
        List of previously accepted half-open query intervals.

    Returns
    -------
    bool
        True if the candidate interval overlaps any accepted interval, otherwise False.
    """
    for accepted_start, accepted_end in accepted_intervals:
        if query_start < accepted_end and accepted_start < query_end:
            return True
    return False


# --------------------------------------------------------------------------------------
# PAF Segmentation
# --------------------------------------------------------------------------------------
def resolve_paf_overlaps(
    paf_file: Path,
    output_file: Path,
) -> Path:
    """
    Resolve overlapping filtered PAF alignments into non-overlapping query segmentation.

    Parameters
    ----------
    paf_file : Path
        Path to the filtered input PAF file.
    output_file : Path
        Path where the resolved non-overlapping PAF file will be written.

    Raises
    ------
    FileNotFoundError
        Raised if the input PAF file does not exist.
    ValueError
        Raised if the PAF file contains malformed rows.

    Returns
    -------
    Path
        Written resolved PAF file path.
    """
    if not paf_file.is_file():
        raise FileNotFoundError(f"PAF file not found: {paf_file}")

    logger.info("Resolving PAF overlaps from %s into %s", paf_file, output_file)

    grouped_alignments = {}
    total_alignments = 0
    discarded_alignments = 0

    with paf_file.open() as file:
        for line_number, line in enumerate(file, start=1):
            if not line.strip():
                continue
            alignment = _parse_paf_alignment(line, line_number)
            grouped_alignments.setdefault(alignment.query_name, []).append(alignment)
            total_alignments += 1

    accepted_lines = []
    kept_alignments = 0

    for query_name, alignments in grouped_alignments.items():
        accepted_intervals = []
        for alignment in sorted(alignments, key=_rank_alignment):
            if _overlaps_existing(alignment.query_start, alignment.query_end, accepted_intervals):
                discarded_alignments += 1
                logger.debug(
                    "Discarding overlapping alignment for query=%s target=%s interval=[%s,%s)",
                    alignment.query_name,
                    alignment.target_name,
                    alignment.query_start,
                    alignment.query_end,
                )
                continue

            accepted_intervals.append((alignment.query_start, alignment.query_end))
            accepted_lines.append("\t".join(alignment.fields) + "\n")
            kept_alignments += 1

        logger.debug(
            "Resolved query=%s to %s non-overlapping alignments",
            query_name,
            len(accepted_intervals),
        )

    if discarded_alignments > 0:
        logger.warning(
            "Discarded %s overlapping PAF alignments while resolving true-atom segmentation",
            discarded_alignments,
        )

    output_file.parent.mkdir(parents=True, exist_ok=True)
    with output_file.open("w") as file:
        file.writelines(accepted_lines)

    logger.info(
        "Resolved PAF overlaps into %s with %s kept alignments out of %s total alignments",
        output_file,
        kept_alignments,
        total_alignments,
    )

    return output_file


# --------------------------------------------------------------------------------------
# GEESE Validation
# --------------------------------------------------------------------------------------
def validate_non_overlapping_geese(geese_file: Path) -> None:
    """
    Validate that a GEESE file contains non-overlapping half-open intervals per genome.

    Parameters
    ----------
    geese_file : Path
        Path to the input GEESE file.

    Raises
    ------
    FileNotFoundError
        Raised if the input GEESE file does not exist.
    ValueError
        Raised if any genome contains overlapping half-open intervals.

    Returns
    -------
    None
    """
    if not geese_file.is_file():
        raise FileNotFoundError(f"GEESE file not found: {geese_file}")

    df = read_geese(geese_file)
    for genome_name, group in df.groupby("name", sort=False):
        sorted_group = group.sort_values(["start", "end"]).reset_index(drop=True)
        previous_end = None
        previous_start = None
        for _, row in sorted_group.iterrows():
            start = int(row["start"])
            end = int(row["end"])
            if previous_end is not None and start < previous_end:
                raise ValueError(
                    "Overlapping true atoms detected for "
                    f"genome '{genome_name}': [{previous_start}, {previous_end}) and [{start}, {end})."
                )
            previous_start = start
            previous_end = end
