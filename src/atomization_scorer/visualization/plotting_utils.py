"""
plotting_utils.py

Shared helper functions for the interactive atomization visualization.

Functions
---------
normalize_output_format : Validate and normalize the interactive output format.
compute_initial_window  : Compute the initial visible genome window for one plot.
get_atoms_for_genome    : Extract validated atom metadata for one genome.
build_class_color_map   : Assign deterministic colors to atom classes.
pair_atoms              : Match true and predicted atoms by class and interval overlap.
sanitize_output_stem    : Convert a genome name into a filesystem-safe HTML stem.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
from __future__ import annotations
from collections import defaultdict
import logging
from typing import Iterable, TypedDict

import pandas as pd

from atomization_scorer.data_processing.utils import sanitize_path_component

# --------------------------------------------------------------------------------------
# Logging
# --------------------------------------------------------------------------------------
log = logging.getLogger(__name__)

# --------------------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------------------
SUPPORTED_OUTPUT_FORMATS = {"html"}
CLASS_PALETTE = [
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
    "#4e79a7",
    "#f28e2b",
    "#59a14f",
    "#e15759",
    "#76b7b2",
    "#edc948",
    "#b07aa1",
    "#ff9da7",
    "#9c755f",
    "#bab0ab",
]


# --------------------------------------------------------------------------------------
# Types
# --------------------------------------------------------------------------------------
class AtomRecord(TypedDict):
    """
    Normalized atom metadata for one genome and one source.

    Attributes
    ----------
    genome_name : str
        Name of the genome this atom belongs to.
    source : str
        Track origin of the atom, either "true" or "predicted".
    class_id : str
        Atom class identifier.
    atom_number : int
        Display atom number used for labeling and identity checks.
    atom_id : str
        Human-readable atom identifier combining class and atom number.
    start : int
        Half-open interval start position on the genome.
    end : int
        Half-open interval end position on the genome.
    length : int
        Length of the atom interval in bases.
    """
    genome_name: str
    source: str
    class_id: str
    atom_number: int
    atom_id: str
    start: int
    end: int
    length: int


# --------------------------------------------------------------------------------------
# Output Format
# --------------------------------------------------------------------------------------
def normalize_output_format(output_format: str) -> str:
    """
    Validate and normalize the interactive output format.

    Parameters
    ----------
    output_format : str
        Requested visualization output format.

    Raises
    ------
    ValueError
        Raised if the requested format is not supported.

    Returns
    -------
    str
        Normalized lowercase output format.
    """
    normalized = output_format.lower()
    if normalized not in SUPPORTED_OUTPUT_FORMATS:
        supported = ", ".join(sorted(SUPPORTED_OUTPUT_FORMATS))
        raise ValueError(
            f"Unsupported output format '{output_format}'. Supported formats: {supported}."
        )
    return normalized


def compute_initial_window(
    genome_length: int,
    target_rows: int,
    min_bases_per_row: int,
    max_bases_per_row: int,
) -> int:
    """
    Compute the initial visible window size for one interactive genome plot.

    Parameters
    ----------
    genome_length : int
        Total genome length in bases.
    target_rows : int
        Desired number of windows across the genome.
    min_bases_per_row : int
        Minimum initial window size in bases.
    max_bases_per_row : int
        Maximum initial window size in bases.

    Returns
    -------
    int
        Initial visible window size in bases.
    """
    window = max(1, genome_length)
    if target_rows > 0:
        window = max(1, -(-genome_length // target_rows))
    window = max(min_bases_per_row, window)
    window = min(max_bases_per_row, window, genome_length)
    return max(1, window)


# --------------------------------------------------------------------------------------
# Atom Extraction
# --------------------------------------------------------------------------------------
def _parse_column_as_int(df: pd.DataFrame, column: str) -> pd.Series | None:
    """
    Parse a DataFrame column as integers, or return None if the column is absent.

    Parameters
    ----------
    df : pd.DataFrame
        Table to read from.
    column : str
        Column name to parse.

    Raises
    ------
    ValueError
        Raised if the column exists but contains non-integer-compatible values.

    Returns
    -------
    pd.Series or None
        Integer Series if the column exists, None otherwise.
    """
    if column not in df.columns:
        return None
    try:
        return pd.to_numeric(df[column], errors="raise").astype(int)
    except (TypeError, ValueError) as error:
        raise ValueError("Atom number column must contain integer-compatible values.") from error


def _parse_atom_number_columns(genome_df: pd.DataFrame) -> pd.Series | None:
    """
    Try to resolve atom numbers from the atom_number or atom column.

    Parameters
    ----------
    genome_df : pd.DataFrame
        Genome-specific atom table.

    Returns
    -------
    pd.Series or None
        Integer Series if a usable column is found, None otherwise.
    """
    result = _parse_column_as_int(genome_df, "atom_number")
    if result is not None:
        return result

    if "atom" in genome_df.columns:
        parsed = pd.to_numeric(genome_df["atom"], errors="coerce")
        if parsed.notna().all():
            return parsed.astype(int)

    return None


def _resolve_display_atom_numbers(genome_df: pd.DataFrame) -> pd.Series:
    """
    Resolve genome-global atom numbers for display and lookup.

    Parameters
    ----------
    genome_df : pd.DataFrame
        Genome-specific atom table.

    Raises
    ------
    ValueError
        Raised if an explicit atom number column is present but not numeric.

    Returns
    -------
    pd.Series
        One-based atom number per row.
    """
    result = _parse_column_as_int(genome_df, "atom_nr")
    if result is not None:
        return result

    result = _parse_atom_number_columns(genome_df)
    if result is not None:
        return result

    return pd.Series(range(1, len(genome_df) + 1), index=genome_df.index, dtype=int)


def get_atoms_for_genome(
    df: pd.DataFrame,
    genome_name: str,
    genome_length: int,
    label: str,
    source: str,
) -> list[AtomRecord]:
    """
    Extract, validate, and normalize atoms for a single genome.

    Parameters
    ----------
    df : pd.DataFrame
        Table containing atom intervals for one or more genomes.
    genome_name : str
        Name of the genome whose atoms are extracted.
    genome_length : int
        Total genome length in bases.
    label : str
        Label used in warning and error messages.
    source : str
        Source label recorded in the output, typically "true" or "predicted".

    Raises
    ------
    ValueError
        Raised if an interval contains a negative coordinate, does not satisfy
        start < end, ends outside the genome length, or has invalid explicit atom numbers.

    Returns
    -------
    list[AtomRecord]
        Sorted list of normalized atom records for the selected genome.
    """
    genome_df = df.loc[df["name"] == genome_name].copy()
    if genome_df.empty:
        return []

    genome_df["start"] = pd.to_numeric(genome_df["start"], errors="raise").astype(int)
    genome_df["end"] = pd.to_numeric(genome_df["end"], errors="raise").astype(int)
    genome_df["class"] = genome_df["class"].astype(str)
    genome_df["atom_number"] = _resolve_display_atom_numbers(genome_df)

    display_df = genome_df.sort_values(["start", "end", "class", "atom_number"], kind="stable")

    previous_end = None
    atoms = []
    for row in display_df.to_dict(orient="records"):
        start = int(row["start"])
        end = int(row["end"])
        atom_number = int(row["atom_number"])
        class_id = str(row["class"])

        if start < 0 or end < 0:
            raise ValueError(f"{label} interval for genome '{genome_name}' contains a negative coordinate.")
        if start >= end:
            raise ValueError(
                f"{label} interval for genome '{genome_name}' must satisfy start < end, got ({start}, {end})."
            )
        if end > genome_length:
            raise ValueError(
                f"{label} interval for genome '{genome_name}' ends outside genome length {genome_length}: ({start}, {end})."
            )
        if previous_end is not None and start < previous_end:
            log.warning(
                "%s intervals for genome '%s' overlap near (%s, %s); visualization will continue.",
                label,
                genome_name,
                start,
                end,
            )
        previous_end = end

        atoms.append(
            AtomRecord(
                genome_name=genome_name,
                source=source,
                class_id=class_id,
                atom_number=atom_number,
                atom_id=f"{class_id}:{atom_number}",
                start=start,
                end=end,
                length=end - start,
            )
        )

    return atoms


# --------------------------------------------------------------------------------------
# Matching and Styling
# --------------------------------------------------------------------------------------
def build_class_color_map(classes: Iterable[str]) -> dict[str, str]:
    """
    Assign deterministic colors to a set of atom classes.

    Parameters
    ----------
    classes : Iterable[str]
        Atom class identifiers.

    Returns
    -------
    dict[str, str]
        Mapping from class identifier to hex color.
    """
    ordered_classes = sorted({str(class_id) for class_id in classes})
    return {
        class_id: CLASS_PALETTE[index % len(CLASS_PALETTE)]
        for index, class_id in enumerate(ordered_classes)
    }


def pair_atoms(
    true_atoms: list[AtomRecord],
    predicted_atoms: list[AtomRecord],
) -> tuple[list[tuple[AtomRecord, AtomRecord]], list[AtomRecord], list[AtomRecord]]:
    """
    Match true and predicted atoms by class and interval overlap.

    Parameters
    ----------
    true_atoms : list[AtomRecord]
        True atoms for one genome.
    predicted_atoms : list[AtomRecord]
        Predicted atoms for one genome.

    Returns
    -------
    tuple[list[tuple[AtomRecord, AtomRecord]], list[AtomRecord], list[AtomRecord]]
        Matched true/predicted pairs, unmatched true atoms, unmatched predicted atoms.
    """
    matched_pairs = []
    matched_true_signatures = set()
    matched_predicted_signatures = set()

    true_by_class = defaultdict(list)
    predicted_by_class = defaultdict(list)

    for atom in true_atoms:
        true_by_class[atom["class_id"]].append(atom)
    for atom in predicted_atoms:
        predicted_by_class[atom["class_id"]].append(atom)

    for class_id in sorted(set(true_by_class) | set(predicted_by_class)):
        class_true_atoms = sorted(
            true_by_class.get(class_id, []),
            key=lambda a: (a["start"], a["end"]),
        )
        class_predicted_atoms = sorted(
            predicted_by_class.get(class_id, []),
            key=lambda a: (a["start"], a["end"]),
        )

        predicted_start_index = 0
        for true_atom in class_true_atoms:
            while (
                predicted_start_index < len(class_predicted_atoms)
                and class_predicted_atoms[predicted_start_index]["end"] <= true_atom["start"]
            ):
                predicted_start_index += 1

            predicted_index = predicted_start_index
            while (
                predicted_index < len(class_predicted_atoms)
                and class_predicted_atoms[predicted_index]["start"] < true_atom["end"]
            ):
                predicted_atom = class_predicted_atoms[predicted_index]
                if true_atom["start"] < predicted_atom["end"] and predicted_atom["start"] < true_atom["end"]:
                    matched_pairs.append((true_atom, predicted_atom))
                    matched_true_signatures.add(
                        (
                            true_atom["class_id"],
                            true_atom["atom_number"],
                            true_atom["start"],
                            true_atom["end"],
                            true_atom["source"],
                        )
                    )
                    matched_predicted_signatures.add(
                        (
                            predicted_atom["class_id"],
                            predicted_atom["atom_number"],
                            predicted_atom["start"],
                            predicted_atom["end"],
                            predicted_atom["source"],
                        )
                    )
                predicted_index += 1

    unmatched_true = [
        atom for atom in true_atoms
        if (
            atom["class_id"],
            atom["atom_number"],
            atom["start"],
            atom["end"],
            atom["source"],
        ) not in matched_true_signatures
    ]
    unmatched_predicted = [
        atom for atom in predicted_atoms
        if (
            atom["class_id"],
            atom["atom_number"],
            atom["start"],
            atom["end"],
            atom["source"],
        ) not in matched_predicted_signatures
    ]

    matched_pairs.sort(key=lambda pair: (pair[0]["start"], pair[0]["end"], pair[0]["class_id"]))
    unmatched_true.sort(key=lambda a: (a["start"], a["end"], a["class_id"]))
    unmatched_predicted.sort(key=lambda a: (a["start"], a["end"], a["class_id"]))
    return matched_pairs, unmatched_true, unmatched_predicted


# --------------------------------------------------------------------------------------
# Output Naming
# --------------------------------------------------------------------------------------
def sanitize_output_stem(value: str) -> str:
    """
    Convert an arbitrary genome identifier into a filesystem-safe HTML stem.

    Parameters
    ----------
    value : str
        Genome name or identifier.

    Returns
    -------
    str
        Filesystem-safe stem.
    """
    return sanitize_path_component(value, fallback_prefix="genome")
