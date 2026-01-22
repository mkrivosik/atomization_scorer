"""
helpers.py

General helper functions for computing metrics and manipulating interval data.

Functions
---------
_compute_and_write_metrics  : Compute precision, recall, and F1-score from TP, FP, and FN
                              and save them to file.
_compute_metrics            : Compute precision, recall, and F1-score from TP, FP, and FN.
_write_metrics_tsv          : Save a DataFrame to a TSV file.
_create_new_row             : Construct a dictionary representing an interval with assigned status.
_interval_overlap           : Compute fraction of overlap between two intervals relative to their union.
"""

# ---------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------
from pathlib import Path

import pandas as pd

# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

def _compute_and_write_metrics(
    tp: dict[int, int],
    fp: dict[int, int],
    fn: dict[int, int],
    output_file: Path,
    per_class: bool = False
) -> float | list[dict[str, float]]:
    """
    Compute precision, recall, and F1-score from TP, FP, FN.

    If `per_class` is False, the function returns the overall F1-score across all classes.
    If `per_class` is True, it returns F1-score, Precision, and Recall for each class.

    Parameters
    ----------
    tp : Dict[int, int]
        Dictionary with the count of True Positives for each atomization class.
    fp : Dict[int, int]
        Dictionary with the count of False Positives for each atomization class.
    fn : Dict[int, int]
        Dictionary with the count of False Negatives for each atomization class.
    output_file : Path
        Path to the output file where the metrics will be saved.
    per_class : bool, optional, default=False
        If True, returns metrics for each class. Otherwise, computes overall metrics.

    Returns
    -------
    float or List[Dict[str, int | float]]
        If per_class is False, returns overall base-level F1-score between 0.0 and 1.0.
        If per_class is True, returns a list of dictionaries, each containing:
            "Class": int -> atomization class,
            "F1-score": float -> F1-score for that class (0.0 to 1.0).
    """
    # -------- overall output --------
    if not per_class:
        total_tp = sum(tp.values())
        total_fp = sum(fp.values())
        total_fn = sum(fn.values())

        precision, recall, f1 = _compute_metrics(total_tp, total_fp, total_fn)

        out_df = pd.DataFrame(
            [[total_tp, total_fp, total_fn, precision, recall, f1]],
            columns=["TP", "FP", "FN", "Precision", "Recall", "F1-score"]
        )
        _write_metrics_tsv(out_df, output_file)
        return f1

    # -------- per-class output --------
    rows = []
    class_f1_metrics = []
    classes = sorted(set(tp) | set(fp) | set(fn))

    for atom_class in classes:
        _tp = tp.get(atom_class, 0)
        _fp = fp.get(atom_class, 0)
        _fn = fn.get(atom_class, 0)

        precision, recall, f1 = _compute_metrics(_tp, _fp, _fn)
        rows.append([atom_class, _tp, _fp, _fn, precision, recall, f1])
        class_f1_metrics.append({"Class": atom_class, "F1-score": f1})

    out_df = pd.DataFrame(
        rows,
        columns=("Class", "TP", "FP", "FN", "Precision", "Recall", "F1-score")
    )
    _write_metrics_tsv(out_df, output_file)
    return class_f1_metrics


def _compute_metrics(tp: int, fp: int, fn: int) -> tuple[float, float, float]:
    """
    Compute precision, recall, and F1-score from True Positives, False Positives, and False Negatives.

    Parameters
    ----------
    tp : int
        Number of True Positives.
    fp : int
        Number of False Positives.
    fn : int
        Number of False Negatives.

    Returns
    -------
    Tuple[float, float, float]
        precision : float
            Precision value between 0.0 and 1.0.
        recall : float
            Recall value between 0.0 and 1.0.
        f1 : float
            F1-score between 0.0 and 1.0.
    """
    precision = tp / (tp + fp) if tp + fp > 0 else 0.0
    recall = tp / (tp + fn) if tp + fn > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall > 0 else 0.0
    return precision, recall, f1


def _write_metrics_tsv(df: pd.DataFrame, output_path: Path) -> None:
    """
    Write metrics DataFrame to TSV file.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing metrics to be written.
    output_path : Path
        Path to the output TSV file.

    Returns
    -------
    None
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep="\t", index=False)


def _create_new_row(
    row: pd.Series,
    start: int,
    end: int,
    status: str
) -> dict[str, str | int]:
    """
    Construct a dictionary representing an interval with assigned status.

    Parameters
    ----------
    row: pd.Series
        Original row from the predicted or true DataFrame.
    start: int
        Start position of the interval.
    end: int
        End position of the interval.
    status: str
        Status of the interval ("TP", "FP", or "FN").

    Returns
    -------
    dict
        Dictionary containing interval information with status.
    """
    return {
        "name": row["name"],
        "atom_nr": row["atom_nr"],
        "class": row["class"],
        "strand": row["strand"],
        "start": start,
        "end": end,
        "status": status,
    }


def _interval_overlap(start1: int, end1: int, start2: int, end2: int) -> float:
    """
    Compute fraction of overlap between two intervals relative to their union.
    Returns 0.0 if intervals do not overlap.

    Parameters
    ----------
    start1 : int
        Start position of first interval.
    end1 : int
        End position of first interval.
    start2 : int
        Start position of second interval.
    end2 : int
        End position of second interval.

    Returns
    -------
    float
        Fraction of overlap (0.0 to 1.0).
    """
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    if overlap_end < overlap_start:
        return 0.0
    overlap = overlap_end - overlap_start + 1
    union = max(end1, end2) - min(start1, start2) + 1
    return overlap / union
