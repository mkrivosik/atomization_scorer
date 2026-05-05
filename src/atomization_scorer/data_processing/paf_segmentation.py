"""
paf_segmentation.py

Resolve filtered PAF alignments into non-overlapping query segmentation.

Functions
---------
_parse_paf_alignment            : Parse one PAF line into a validated alignment record.
_rank_alignment                 : Build a sorting key for prioritizing candidate alignments.
_compute_surviving_intervals    : Subtract accepted intervals from a candidate interval and
                                  return the remaining sub-intervals.
resolve_paf_overlaps            : Resolve overlapping filtered PAF alignments into valid
                                  non-overlapping query segmentation, trimming candidates
                                  at overlap boundaries before applying a minimum length threshold.
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
from typing import cast

from .geese_reader import read_geese

# --------------------------------------------------------------------------------------
# Logging
# --------------------------------------------------------------------------------------
log = logging.getLogger(__name__)


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
# Interval Subtraction
# --------------------------------------------------------------------------------------
def _compute_surviving_intervals(
    query_start: int,
    query_end: int,
    accepted_intervals: list[tuple[int, int]],
) -> list[tuple[int, int]]:
    """
    Return the portions of [query_start, query_end) not covered by any accepted interval.

    Parameters
    ----------
    query_start : int
        Start position of the candidate half-open interval.
    query_end : int
        End position of the candidate half-open interval.
    accepted_intervals : list[tuple[int, int]]
        Already-accepted half-open intervals to subtract from the candidate.

    Returns
    -------
    list[tuple[int, int]]
        Remaining sub-intervals after removing all accepted coverage. Empty when the
        candidate is fully contained within accepted intervals.
    """
    free_regions = [(query_start, query_end)]
    for accepted_start, accepted_end in accepted_intervals:
        if accepted_start >= query_end or accepted_end <= query_start:
            continue
        next_free = []
        for free_start, free_end in free_regions:
            if accepted_start >= free_end or accepted_end <= free_start:
                next_free.append((free_start, free_end))
            else:
                if free_start < accepted_start:
                    next_free.append((free_start, accepted_start))
                if accepted_end < free_end:
                    next_free.append((accepted_end, free_end))
        free_regions = next_free
    return free_regions


# --------------------------------------------------------------------------------------
# PAF Segmentation
# --------------------------------------------------------------------------------------
def resolve_paf_overlaps(
    paf_file: Path,
    output_file: Path,
    minimum_alignment_length: int = 1,
) -> Path:
    """
    Resolve overlapping filtered PAF alignments into non-overlapping query segmentation.

    Non-primary alignments (any alignment_type other than "P") are discarded before
    overlap resolution begins. Only primary alignments compete for query-sequence space.

    Parameters
    ----------
    paf_file : Path
        Path to the filtered input PAF file.
    output_file : Path
        Path where the resolved non-overlapping PAF file will be written.
    minimum_alignment_length : int, optional, default=1
        Minimum length in bases that a trimmed sub-interval must reach to be kept.

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

    Notes
    -----
    Trimmed rows only update fields[2] (query_start) and fields[3] (query_end).
    Alignment length, match count, and target coordinates retain the original values
    and no longer describe the surviving sub-interval (because that information will
    not be needed downstream).
    """
    if not paf_file.is_file():
        raise FileNotFoundError(f"PAF file not found: {paf_file}")

    log.info("Resolving PAF overlaps from %s into %s", paf_file, output_file)

    grouped_alignments = {}
    total_alignments = 0
    discarded_alignments = 0
    non_primary_count = 0

    with paf_file.open() as file:
        for line_number, line in enumerate(file, start=1):
            if not line.strip():
                continue
            alignment = _parse_paf_alignment(line, line_number)
            if alignment.alignment_type != "P":
                non_primary_count += 1
                continue
            grouped_alignments.setdefault(alignment.query_name, []).append(alignment)
            total_alignments += 1

    if non_primary_count > 0:
        log.info("Discarded %s non-primary alignments before overlap resolution", non_primary_count)

    accepted_lines = []
    kept_alignments = 0

    for query_name, alignments in grouped_alignments.items():
        accepted_intervals = []
        for alignment in sorted(alignments, key=_rank_alignment):
            surviving = _compute_surviving_intervals(
                alignment.query_start,
                alignment.query_end,
                accepted_intervals,
            )
            kept_any = False
            for trim_start, trim_end in surviving:
                if trim_end - trim_start >= minimum_alignment_length:
                    trimmed_fields = list(alignment.fields)
                    trimmed_fields[2] = str(trim_start)
                    trimmed_fields[3] = str(trim_end)
                    accepted_intervals.append((trim_start, trim_end))
                    accepted_lines.append("\t".join(trimmed_fields) + "\n")
                    kept_alignments += 1
                    kept_any = True
                else:
                    log.debug(
                        "Discarding trimmed interval [%s, %s) for query=%s target=%s: "
                        "length %s < minimum_alignment_length %s",
                        trim_start,
                        trim_end,
                        alignment.query_name,
                        alignment.target_name,
                        trim_end - trim_start,
                        minimum_alignment_length,
                    )
            if not kept_any:
                discarded_alignments += 1
                log.debug(
                    "Fully discarding alignment for query=%s target=%s interval=[%s,%s): "
                    "no surviving portion >= minimum_alignment_length=%s",
                    alignment.query_name,
                    alignment.target_name,
                    alignment.query_start,
                    alignment.query_end,
                    minimum_alignment_length,
                )

        log.debug(
            "Resolved query=%s to %s non-overlapping alignments",
            query_name,
            len(accepted_intervals),
        )

    if discarded_alignments > 0:
        log.warning(
            "Fully discarded %s PAF alignments during overlap resolution "
            "(no surviving portion >= minimum_alignment_length=%s)",
            discarded_alignments,
            minimum_alignment_length,
        )

    output_file.parent.mkdir(parents=True, exist_ok=True)
    with output_file.open("w") as file:
        file.writelines(accepted_lines)

    log.info(
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
    for genome_name_key, group in df.groupby("name", sort=False):
        genome_name = cast(str, genome_name_key)
        sorted_group = group.sort_values(["start", "end"]).reset_index(drop=True)
        starts = [int(value) for value in sorted_group["start"]]
        ends = [int(value) for value in sorted_group["end"]]
        for (previous_start, previous_end), (start, end) in zip(zip(starts, ends), zip(starts[1:], ends[1:])):
            if start < previous_end:
                raise ValueError(
                    "Overlapping true atoms detected for "
                    f"genome '{genome_name}': [{previous_start}, {previous_end}) and [{start}, {end})."
                )
