"""
interval_metrics.py

Compute interval-level alignment metrics (TP, FP, FN, Precision, Recall, F1-score) for genome atomization.

Functions
---------
compute_interval_level_metrics  : Computes per-interval alignment metrics and F1-score.
_scan_intervals_interval_level  : Scan predicted and true intervals to count TP, FP, and FN at interval-level.
_write_interval_status          : Write predicted and true intervals with status (TP/FP/FN) to TSV file.
"""

# ---------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------
from pathlib import Path

import pandas as pd
from intervaltree import IntervalTree

from atomization_scorer.data_processing import read_geese

from .helpers import _compute_and_write_metrics, _interval_overlap, _write_metrics_tsv

# ---------------------------------------------------------------------
# Interval-Level Metrics
# ---------------------------------------------------------------------

def compute_interval_level_metrics(
    predicted_geese: Path,
    true_geese: Path,
    output_directory: Path,
    per_class: bool = False,
    min_overlap_ratio: float = 0.8
) -> float | list[dict[str, int | float]]:
    """
    Compute interval-level True Positives, False Positives, False Negatives,
    Precision, Recall, and F1-score.

    Parameters
    ----------
    predicted_geese : Path
        Input GEESE file containing the predicted atomization.
    true_geese : Path
        Input GEESE file containing the true (gold standard) atomization.
    output_directory : Path
       Path to the output directory where results are stored.
    per_class : bool, optional, default=False
        If True, compute metrics per atomization class; otherwise overall.
    min_overlap_ratio : float, optional, default=0.8
        Minimum overlap ratio for interval-level scoring.

    Raises
    ------
    FileNotFoundError
        Raised if predicted_geese or true_geese do not exist.

    Returns
    -------
    float or List[Dict[str, int | float]]
        If per_class is False, returns overall interval-level F1-score between 0.0 and 1.0.
        If per_class is True, returns a list of dictionaries, each containing:
            "Class": int -> atomization class,
            "F1-score": float -> F1-score for that class (0.0 to 1.0).
    """
    if not predicted_geese.is_file():
        raise FileNotFoundError(f"Predicted GEESE file not found: {predicted_geese}")
    if not true_geese.is_file():
        raise FileNotFoundError(f"True GEESE file not found: {true_geese}")

    output_directory.mkdir(parents=True, exist_ok=True)
    if per_class:
        output_file = output_directory / "interval_metrics_per_class.tsv"
    else:
        output_file = output_directory / "interval_metrics_overall.tsv"

    predicted_df = read_geese(geese_file=predicted_geese)
    true_df = read_geese(geese_file=true_geese)

    predicted_df = predicted_df.sort_values(["start", "end"]).reset_index(drop=True)
    true_df = true_df.sort_values(["start", "end"]).reset_index(drop=True)

    tp, fp, fn = _scan_intervals_interval_level(
        predicted_df=predicted_df,
        true_df=true_df,
        output_directory=output_directory,
        min_overlap_ratio=min_overlap_ratio
    )

    return _compute_and_write_metrics(
        tp=tp, fp=fp, fn=fn,
        output_file=output_file,
        per_class=per_class
    )


def _scan_intervals_interval_level(
    predicted_df: pd.DataFrame,
    true_df: pd.DataFrame,
    output_directory: Path,
    min_overlap_ratio: float = 0.8
) -> tuple[
    dict[int, int],
    dict[int, int],
    dict[int, int]
]:
    """
    Scan predicted and true intervals to compute True Positives (TP), False Positives (FP), and False Negatives (FN)
    per atomization class. Counts are computed at interval-level.

    Parameters
    ----------
    predicted_df : pd.DataFrame
        DataFrame containing predicted intervals.
    true_df : pd.DataFrame
        DataFrame containing true intervals.
    output_directory : Path
        Path to the output directory where results are stored.
    min_overlap_ratio : float, optional, default=0.8
        Minimum overlap ratio for counting a predicted interval as True Positive.

    Returns
    -------
    Tuple[Dict[int, int], Dict[int, int], Dict[int, int]]
        tp : Dict[int, int]
            True Positives per atom class.
        fp : Dict[int, int]
            False Positives per atom class.
        fn : Dict[int, int]
            False Negatives per atom class.
    """
    tp = {}  # True Positives
    fp = {}  # False Positives
    fn = {}  # False Negatives
    tp_indexes = set()

    all_sequences = set(predicted_df["name"]).union(set(true_df["name"]))

    for sequence in all_sequences:
        predicted_sequence = predicted_df[predicted_df["name"] == sequence].reset_index(drop=True)
        true_sequence = true_df[true_df["name"] == sequence].reset_index(drop=True)

        tree = IntervalTree()
        for index, row in true_sequence.iterrows():
            tree[row["start"]:row["end"] + 1] = (index, int(row["class"]))

        true_used = [False] * len(true_sequence)

        for predicted_index, predicted in predicted_sequence.iterrows():
            overlaps = tree[predicted["start"]:predicted["end"] + 1]
            matched = False

            for interval in overlaps:
                true_index, true_class = interval.data
                if not true_used[true_index] and int(predicted["class"]) == true_class:
                    overlap_ratio = _interval_overlap(
                        predicted["start"], predicted["end"],
                        true_sequence.iloc[true_index]["start"], true_sequence.iloc[true_index]["end"]
                    )
                    if overlap_ratio >= min_overlap_ratio:
                        atom_class = int(predicted["class"])
                        tp[atom_class] = tp.get(atom_class, 0) + 1
                        true_used[true_index] = True
                        matched = True
                        tp_indexes.add(predicted_index)
                        break

            if not matched:
                atom_class = int(predicted["class"])
                fp[atom_class] = fp.get(atom_class, 0) + 1

        for index, used in enumerate(true_used):
            if not used:
                atom_class = int(true_sequence.iloc[index]["class"])
                fn[atom_class] = fn.get(atom_class, 0) + 1

    predicted_file = output_directory / "interval_predicted_status.tsv"
    true_file = output_directory / "interval_true_status.tsv"
    _write_interval_status(
        predicted_df=predicted_df,
        true_df=true_df,
        tp_indexes=tp_indexes,
        predicted_output_path=predicted_file,
        true_output_path=true_file
    )

    return tp, fp, fn


def _write_interval_status(
    predicted_df: pd.DataFrame,
    true_df: pd.DataFrame,
    tp_indexes: set[int],
    predicted_output_path: Path,
    true_output_path: Path
) -> None:
    """
    Write predicted and true intervals to TSV files, adding a "status" column to indicate whether each
    interval is a True Positive (TP), False Positive (FP), or False Negative (FN).

    Parameters
    ----------
    predicted_df : pd.DataFrame
        DataFrame containing the predicted intervals.
    true_df : pd.DataFrame
        DataFrame containing the true intervals.
    tp_indexes : Set[int]
        Set of indexes of predicted intervals that were classified as True Positives;
        same indexes are used to mark True Positives in the true intervals as well.
    predicted_output_path : Path
        Path where the predicted intervals TSV file will be saved.
    true_output_path : Path
        Path where the true intervals TSV file will be saved.

    Returns
    -------
    None
    """
    # -------- predicted intervals --------
    predicted_out = predicted_df.copy()
    predicted_status = []

    for i in range(len(predicted_df)):
        if i in tp_indexes:
            predicted_status.append("TP")
        else:
            predicted_status.append("FP")

    predicted_out["status"] = predicted_status
    _write_metrics_tsv(df=predicted_out, output_path=predicted_output_path)

    # -------- true intervals --------
    true_out = true_df.copy()
    true_status = []

    for i in range(len(true_df)):
        if i in tp_indexes:
            true_status.append("TP")
        else:
            true_status.append("FN")

    true_out["status"] = true_status
    _write_metrics_tsv(df=true_out, output_path=true_output_path)
