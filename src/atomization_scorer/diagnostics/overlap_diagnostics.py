"""
overlap_diagnostics.py

Provides utilities for diagnosing overlapping filtered PAF alignments before
true-atom overlap resolution.

Functions
---------
_parse_paf_alignment   : Parse one PAF row into a structured diagnostic alignment record.
_classify_overlap      : Classify a query-side overlap between two diagnostic alignments.
_detect_overlaps       : Detect all overlapping alignment pairs within one query group.
_write_tsv             : Write one tab-delimited diagnostics table.
_write_overlaps_json   : Write reported overlaps into a nested JSON file.
diagnose_paf_overlaps  : Generate overlap reports and dotplot FASTA inputs from filtered PAF.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
from __future__ import annotations
import json
from collections import defaultdict
from dataclasses import dataclass
import logging
from pathlib import Path

# --------------------------------------------------------------------------------------
# Logging
# --------------------------------------------------------------------------------------
log = logging.getLogger(__name__)

# --------------------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------------------
OVERLAP_FIELDNAMES = [
    "query_name",
    "anchor_atom",
    "partner_atom",
    "anchor_query_start",
    "anchor_query_end",
    "partner_query_start",
    "partner_query_end",
    "anchor_atom_length",
    "partner_atom_length",
    "anchor_target_start",
    "anchor_target_end",
    "partner_target_start",
    "partner_target_end",
    "anchor_strand",
    "partner_strand",
    "overlap_start",
    "overlap_end",
    "overlap_length",
    "overlap_class",
]


# --------------------------------------------------------------------------------------
# Diagnostic Alignment Record
# --------------------------------------------------------------------------------------
@dataclass(frozen=True)
class DiagnosticAlignment:
    """
    Store parsed PAF alignment fields required for overlap diagnostics.

    Attributes
    ----------
    query_name : str
        Name of the query genome sequence.
    query_length : int
        Total length of the query genome sequence.
    query_start : int
        Start position of the aligned interval on the query genome.
    query_end : int
        End position of the aligned interval on the query genome.
    strand : str
        Strand orientation of the alignment ("+" or "-").
    target_name : str
        Name of the target representative atom sequence.
    target_length : int
        Total length of the target representative atom.
    target_start : int
        Start position of the aligned interval on the target atom.
    target_end : int
        End position of the aligned interval on the target atom.
    matches : int
        Number of matching bases in the alignment.
    alignment_length : int
        Length of the aligned block.
    mapq : int
        Minimap2 mapping quality score.
    similarity : float
        Estimated alignment similarity.
    """
    query_name: str
    query_length: int
    query_start: int
    query_end: int
    strand: str
    target_name: str
    target_length: int
    target_start: int
    target_end: int
    matches: int
    alignment_length: int
    mapq: int
    similarity: float


# --------------------------------------------------------------------------------------
# PAF Alignment Parsing
# --------------------------------------------------------------------------------------
def _parse_paf_alignment(line: str, line_number: int) -> DiagnosticAlignment:
    """
    Parse a PAF row into a structured alignment record for overlap diagnostics.

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
    DiagnosticAlignment
        Parsed diagnostic alignment record with query and target coordinates.
    """
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 12:
        raise ValueError(
            f"Malformed PAF row at line {line_number}: expected at least 12 tab-separated fields."
        )

    try:
        query_length = int(fields[1])
        query_start = int(fields[2])
        query_end = int(fields[3])
        target_length = int(fields[6])
        target_start = int(fields[7])
        target_end = int(fields[8])
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
    if target_start >= target_end:
        raise ValueError(
            f"Malformed PAF row at line {line_number}: target start must be smaller than target end."
        )

    similarity = matches / alignment_length if alignment_length > 0 else 0.0
    for field in fields[12:]:
        if field.startswith("dv:f:"):
            similarity = 1.0 - float(field.split(":")[2])
            break

    return DiagnosticAlignment(
        query_name=fields[0],
        query_length=query_length,
        query_start=query_start,
        query_end=query_end,
        strand=fields[4],
        target_name=fields[5],
        target_length=target_length,
        target_start=target_start,
        target_end=target_end,
        matches=matches,
        alignment_length=alignment_length,
        mapq=mapq,
        similarity=similarity,
    )


# --------------------------------------------------------------------------------------
# Overlap Classification
# --------------------------------------------------------------------------------------
def _classify_overlap(
    overlap_start: int,
    overlap_end: int,
    left: DiagnosticAlignment,
    right: DiagnosticAlignment,
) -> str:
    """
    Classify the query-side overlap between two diagnostic alignments.

    Parameters
    ----------
    overlap_start : int
        Start position of the overlap on the query genome.
    overlap_end : int
        End position of the overlap on the query genome.
    left : DiagnosticAlignment
        First diagnostic alignment in the overlapping pair.
    right : DiagnosticAlignment
        Second diagnostic alignment in the overlapping pair.

    Returns
    -------
    str
        Overlap class label:
        "both_edge" if the overlap touches an edge of both alignments,
        "mixed_edge_internal" otherwise.
    """
    left_edge = overlap_start == left.query_start or overlap_end == left.query_end
    right_edge = overlap_start == right.query_start or overlap_end == right.query_end
    if left_edge and right_edge:
        return "both_edge"
    return "mixed_edge_internal"


# --------------------------------------------------------------------------------------
# Overlap Detection
# --------------------------------------------------------------------------------------
def _detect_overlaps(alignments: list[DiagnosticAlignment]) -> list[dict[str, str]]:
    """
    Detect all overlapping alignment pairs within one query group.

    Parameters
    ----------
    alignments : list[DiagnosticAlignment]
        Diagnostic alignments belonging to the same query genome.

    Returns
    -------
    list[dict[str, str]]
        List of overlap records containing atom IDs, query coordinates,
        target coordinates, strand information, overlap coordinates,
        overlap length, and overlap class.
    """
    overlaps = []
    sorted_alignments = sorted(
        alignments,
        key=lambda alignment: (
            alignment.query_start,
            alignment.query_end,
            alignment.target_name,
        ),
    )

    for left_index, left in enumerate(sorted_alignments):
        for right in sorted_alignments[left_index + 1:]:
            if right.query_start >= left.query_end:
                break
            if left.target_name == right.target_name:
                continue

            overlap_start = max(left.query_start, right.query_start)
            overlap_end = min(left.query_end, right.query_end)
            if overlap_start >= overlap_end:
                continue

            anchor, partner = sorted([left, right], key=lambda item: item.target_name)
            overlaps.append(
                {
                    "query_name": left.query_name,
                    "anchor_atom": anchor.target_name,
                    "partner_atom": partner.target_name,
                    "anchor_query_start": str(anchor.query_start),
                    "anchor_query_end": str(anchor.query_end),
                    "partner_query_start": str(partner.query_start),
                    "partner_query_end": str(partner.query_end),
                    "anchor_atom_length": str(anchor.target_length),
                    "partner_atom_length": str(partner.target_length),
                    "anchor_target_start": str(anchor.target_start),
                    "anchor_target_end": str(anchor.target_end),
                    "partner_target_start": str(partner.target_start),
                    "partner_target_end": str(partner.target_end),
                    "anchor_strand": anchor.strand,
                    "partner_strand": partner.strand,
                    "overlap_start": str(overlap_start),
                    "overlap_end": str(overlap_end),
                    "overlap_length": str(overlap_end - overlap_start),
                    "overlap_class": _classify_overlap(overlap_start, overlap_end, anchor, partner),
                }
            )

    return overlaps


# --------------------------------------------------------------------------------------
# Table/JSON Writers
# --------------------------------------------------------------------------------------
def _write_tsv(
    output_file: Path,
    rows: list[dict[str, str]],
    fieldnames: list[str],
    *,
    pad_columns: bool = True,
    header_separator: bool = False,
    header_labels: dict[str, str] | None = None,
) -> None:
    """
    Write a visually padded tab-delimited table.

    Parameters
    ----------
    output_file : Path
        Path where the TSV file should be written.
    rows : list[dict[str, str]]
        List of row dictionaries to be written into the TSV file.
    fieldnames : list[str]
        Ordered list of column names used for the TSV header and row layout.
    pad_columns : bool, optional, default=True
        Whether to pad columns to a fixed width for easier terminal viewing.
    header_separator : bool, optional, default=False
        Whether to insert one dashed separator row after the header.
    header_labels : dict[str, str] or None, optional, default=None
        Optional display labels to use in the written header row.

    Returns
    -------
    None
    """
    header_names = [
        header_labels.get(fieldname, fieldname)
        if header_labels else fieldname
        for fieldname in fieldnames
    ]
    column_widths = {
        fieldname: max(
            len(header_names[index]),
            max((len(str(row.get(fieldname, ""))) for row in rows), default=0),
        )
        for index, fieldname in enumerate(fieldnames)
    }

    output_file.parent.mkdir(parents=True, exist_ok=True)
    with output_file.open("w", newline="") as file:
        if pad_columns:
            header_line = "\t".join(
                header_names[index].ljust(column_widths[fieldname])
                for index, fieldname in enumerate(fieldnames)
            )
        else:
            header_line = "\t".join(header_names)
        file.write(header_line + "\n")
        if header_separator:
            if pad_columns:
                separator_line = "\t".join("-" * column_widths[fieldname] for fieldname in fieldnames)
            else:
                separator_line = "\t".join("-" * len(header_names[index]) for index, _fieldname in enumerate(fieldnames))
            file.write(separator_line + "\n")
        for row in rows:
            if pad_columns:
                line = "\t".join(
                    str(row.get(fieldname, "")).ljust(column_widths[fieldname])
                    for fieldname in fieldnames
                )
            else:
                line = "\t".join(str(row.get(fieldname, "")) for fieldname in fieldnames)
            file.write(line + "\n")


def _write_overlaps_json(output_file: Path, overlaps: list[dict[str, str]]) -> None:
    """
    Write reported overlaps into a nested JSON file.

    Parameters
    ----------
    output_file : Path
        Path where the nested overlap JSON file should be written.
    overlaps : list[dict[str, str]]
        Flat overlap dictionaries produced by the diagnostics workflow.

    Returns
    -------
    None
    """
    nested_overlaps = [
        {
            "query_name": overlap["query_name"],
            "anchor": {
                "atom": overlap["anchor_atom"],
                "representative_atom_length": int(overlap["anchor_atom_length"]),
                "query_start": int(overlap["anchor_query_start"]),
                "query_end": int(overlap["anchor_query_end"]),
                "target_start": int(overlap["anchor_target_start"]),
                "target_end": int(overlap["anchor_target_end"]),
                "strand": overlap["anchor_strand"],
            },
            "partner": {
                "atom": overlap["partner_atom"],
                "representative_atom_length": int(overlap["partner_atom_length"]),
                "query_start": int(overlap["partner_query_start"]),
                "query_end": int(overlap["partner_query_end"]),
                "target_start": int(overlap["partner_target_start"]),
                "target_end": int(overlap["partner_target_end"]),
                "strand": overlap["partner_strand"],
            },
            "overlap": {
                "start": int(overlap["overlap_start"]),
                "end": int(overlap["overlap_end"]),
                "length": int(overlap["overlap_length"]),
                "class": overlap["overlap_class"],
            },
        }
        for overlap in overlaps
    ]

    output_file.parent.mkdir(parents=True, exist_ok=True)
    with output_file.open("w") as file:
        json.dump(nested_overlaps, file, indent=2)
        file.write("\n")


# --------------------------------------------------------------------------------------
# Overlap Diagnostics
# --------------------------------------------------------------------------------------
def diagnose_paf_overlaps(
    paf_file: Path,
    representatives_fasta: Path,
    output_directory: Path,
    minimum_report_overlap_length: int = 0,
    minimum_plot_overlap_length: int = 0,
    include_reverse: bool = False,
    run_dotter: bool = True,
    dotter_extra_args: list[str] | None = None,
) -> Path:
    """
    Generate overlap reports and anchor FASTA files from a filtered PAF file.

    Parameters
    ----------
    paf_file : Path
        Filtered PAF file inspected before overlap resolution.
    representatives_fasta : Path
        FASTA file containing representative atom sequences.
    output_directory : Path
        Directory where diagnostic outputs should be written.
    minimum_report_overlap_length : int, optional, default=0
        Overlaps at or below this threshold are counted in the summary but excluded from
        anchor-level reports.
    minimum_plot_overlap_length : int, optional, default=0
        Overlaps above this threshold get FASTA files for dotplot inspection.
    include_reverse : bool, optional, default=False
        Whether to duplicate anchor-partner diagnostics from both anchor perspectives.
    run_dotter : bool, optional, default=True
        Whether to run Dotter immediately after generating anchor FASTA inputs.
    dotter_extra_args : list[str] or None, optional, default=None
        Additional command-line arguments passed to Dotter when run_dotter is enabled.

    Raises
    ------
    FileNotFoundError
        Raised if the filtered PAF file or representative FASTA file does not exist.

    Returns
    -------
    Path
        Diagnostic output directory.
    """
    if not paf_file.is_file():
        raise FileNotFoundError(f"PAF file not found: {paf_file}")
    if not representatives_fasta.is_file():
        raise FileNotFoundError(f"Representative FASTA file not found: {representatives_fasta}")

    output_directory.mkdir(parents=True, exist_ok=True)
    anchors_directory = output_directory / "anchors"
    anchors_directory.mkdir(parents=True, exist_ok=True)

    grouped_alignments = defaultdict(list)
    total_alignments = 0

    with paf_file.open() as file:
        for line_number, line in enumerate(file, start=1):
            if not line.strip():
                continue
            alignment = _parse_paf_alignment(line, line_number)
            grouped_alignments[alignment.query_name].append(alignment)
            total_alignments += 1

    all_overlaps = []
    for alignments in grouped_alignments.values():
        all_overlaps.extend(_detect_overlaps(alignments))

    reported_overlaps = [
        overlap
        for overlap in all_overlaps
        if int(overlap["overlap_length"]) > minimum_report_overlap_length
    ]

    _write_overlaps_json(
        output_file=output_directory / "overlaps.json",
        overlaps=reported_overlaps,
    )

    anchor_genome_stats = defaultdict(
        lambda: {
            "partners": set(),
            "count": 0,
            "max_overlap_length": 0,
            "both_edge": 0,
            "mixed_edge_internal": 0,
        }
    )
    for overlap in reported_overlaps:
        overlap_length = int(overlap["overlap_length"])
        overlap_class = overlap["overlap_class"]
        anchor_partner_pairs = [(overlap["anchor_atom"], overlap["partner_atom"])]
        if include_reverse:
            anchor_partner_pairs.append((overlap["partner_atom"], overlap["anchor_atom"]))
        for anchor_atom, partner_atom in anchor_partner_pairs:
            stats = anchor_genome_stats[(overlap["query_name"], anchor_atom)]
            stats["partners"].add(partner_atom)
            stats["count"] += 1
            stats["max_overlap_length"] = max(int(stats["max_overlap_length"]), overlap_length)
            stats[overlap_class] += 1

    anchor_genome_summary_rows = [
        {
            "query_name": query_name,
            "anchor_atom": anchor_atom,
            "n_partners": str(len(stats["partners"])),
            "n_overlaps": str(int(stats["count"])),
            "max_overlap_length": str(int(stats["max_overlap_length"])),
            "n_both_edge": str(int(stats["both_edge"])),
            "n_one_eaten": str(int(stats["mixed_edge_internal"])),
        }
        for (query_name, anchor_atom), stats in sorted(anchor_genome_stats.items())
    ]
    _write_tsv(
        output_file=output_directory / "anchor_genome_summary.tsv",
        rows=anchor_genome_summary_rows,
        fieldnames=[
            "query_name",
            "anchor_atom",
            "n_partners",
            "n_overlaps",
            "max_overlap_length",
            "n_both_edge",
            "n_one_eaten",
        ],
        header_separator=True,
    )

    summary_rows = [
        {"metric": "total_filtered_alignments", "value": str(total_alignments)},
        {"metric": "total_overlapping_pairs", "value": str(len(all_overlaps))},
        {
            "metric": "pairs_overlap_1_to_9bp",
            "value": str(
                sum(
                    1
                    for overlap in all_overlaps
                    if 0 < int(overlap["overlap_length"]) < 10
                )
            ),
        },
        {
            "metric": "pairs_overlap_10_to_99bp",
            "value": str(
                sum(
                    1
                    for overlap in all_overlaps
                    if 10 <= int(overlap["overlap_length"]) < 100
                )
            ),
        },
        {
            "metric": "pairs_overlap_ge_100bp",
            "value": str(
                sum(
                    1
                    for overlap in all_overlaps
                    if 100 <= int(overlap["overlap_length"]) < 1000
                )
            ),
        },
        {
            "metric": "pairs_overlap_1000_to_9999bp",
            "value": str(
                sum(
                    1
                    for overlap in all_overlaps
                    if 1000 <= int(overlap["overlap_length"]) < 10000
                )
            ),
        },
        {
            "metric": "pairs_overlap_ge_10000bp",
            "value": str(
                sum(
                    1
                    for overlap in all_overlaps
                    if 10000 <= int(overlap["overlap_length"])
                )
            ),
        },
        {
            "metric": "pairs_both_edge",
            "value": str(
                sum(
                    1
                    for overlap in all_overlaps
                    if overlap["overlap_class"] == "both_edge"
                )
            ),
        },
        {
            "metric": "pairs_one_eaten",
            "value": str(
                sum(
                    1
                    for overlap in all_overlaps
                    if overlap["overlap_class"] == "mixed_edge_internal"
                )
            ),
        },
        {
            "metric": "anchors_with_reported_overlaps",
            "value": str(len({row["anchor_atom"] for row in anchor_genome_summary_rows})),
        },
        {
            "metric": "query_anchor_pairs_with_reported_overlaps",
            "value": str(len(anchor_genome_summary_rows)),
        },
        {
            "metric": "max_overlap_length",
            "value": str(
                max((int(overlap["overlap_length"]) for overlap in all_overlaps), default=0)
            ),
        },
    ]
    _write_tsv(
        output_file=output_directory / "summary.tsv",
        rows=summary_rows,
        fieldnames=["metric", "value"],
        header_separator=True,
    )

    from .dotplot_inputs import write_anchor_dotplot_fastas

    write_anchor_dotplot_fastas(
        overlaps=reported_overlaps,
        representatives_fasta=representatives_fasta,
        output_directory=anchors_directory,
        minimum_overlap_length=minimum_plot_overlap_length,
        include_reverse=include_reverse,
    )

    if run_dotter:
        from .dotter_runner import run_dotter_for_anchors

        log.info(
            "Running Dotter for anchor FASTA inputs in %s",
            anchors_directory,
        )
        run_dotter_for_anchors(
            anchors_directory=anchors_directory,
            extra_args=dotter_extra_args,
        )

    log.info(
        "Generated overlap diagnostics at %s with %s reported overlaps",
        output_directory,
        len(reported_overlaps),
    )
    return output_directory
