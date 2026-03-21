"""
base_metrics.py

Compute base-level alignment metrics (TP, FP, FN, Precision, Recall, F1-score) for genome atomization.

Functions
---------
compute_base_level_metrics  : Computes per-base alignment metrics and F1-score.
_scan_intervals_base_level  : Scan predicted and true intervals to count TP, FP, and FN at base-level.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
from pathlib import Path

import pandas as pd

from atomization_scorer.data_processing import read_geese

from .helpers import _compute_and_write_metrics, _create_new_row, _write_metrics_tsv

# --------------------------------------------------------------------------------------
# Base-Level Metrics
# --------------------------------------------------------------------------------------
EVENT_PRIORITY = {
    "predicted_end": 0,
    "true_end": 1,
    "predicted_start": 2,
    "true_start": 3,
}


def compute_base_level_metrics(
    predicted_geese: Path,
    true_geese: Path,
    output_directory: Path,
    per_class: bool = False
) -> float | list[dict[str, int | float]]:
    """
    Compute base-level True Positives, False Positives, False Negatives,
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

    Raises
    ------
    FileNotFoundError
        Raised if predicted_geese or true_geese do not exist.

    Returns
    -------
    float or List[Dict[str, int | float]]
        If per_class is False, returns overall base-level F1-score between 0.0 and 1.0.
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
        output_file = output_directory / "base_metrics_per_class.tsv"
    else:
        output_file = output_directory / "base_metrics_overall.tsv"

    predicted_df = read_geese(geese_file=predicted_geese)
    true_df = read_geese(geese_file=true_geese)

    tp, fp, fn = _scan_intervals_base_level(
        predicted_df=predicted_df,
        true_df=true_df,
        output_directory=output_directory,
    )

    return _compute_and_write_metrics(
        tp=tp, fp=fp, fn=fn,
        output_file=output_file,
        per_class=per_class
    )


def _scan_intervals_base_level(
        predicted_df: pd.DataFrame,
        true_df: pd.DataFrame,
        output_directory: Path
) -> tuple[
    dict[int, int],
    dict[int, int],
    dict[int, int]
]:
    """
    Scan predicted and true intervals to compute True Positives (TP), False Positives (FP), and False Negatives (FN)
    per atomization class. Counts are computed at base-level.

    Algorithm
    ---------
    Sweep Line Algorithm (from computational geometry) is used to efficiently calculate per-base True Positives (TP),
    False Positives (FP), and False Negatives (FN) by tracking active predicted and true intervals along the sequence.

    Parameters
    ----------
    predicted_df: pd.DataFrame
        Predicted atomization intervals.
    true_df: pd.DataFrame
        True (gold standard) atomization intervals.
    output_directory: Path
        Path to the output directory where results are stored.

    Returns
    -------
    Tuple[Dict[int, int], Dict[int, int], Dict[int, int]]
        tp : Dict[int, int]
            Number of True Positive bases per Atomization class.
        fp : Dict[int, int]
            Number of False Positive bases per Atomization class.
        fn : Dict[int, int]
            Number of False Negative bases per Atomization class.
    """
    tp = {}  # True Positives
    fp = {}  # False Positives
    fn = {}  # False Negatives
    predicted_intervals = []
    true_intervals = []
    predicted_groups = {name: group for name, group in predicted_df.groupby("name", sort=False)}
    true_groups = {name: group for name, group in true_df.groupby("name", sort=False)}
    all_sequences = set(predicted_groups).union(set(true_groups))

    for sequence in all_sequences:
        predicted_sequence = predicted_groups.get(sequence, predicted_df.iloc[0:0])
        true_sequence = true_groups.get(sequence, true_df.iloc[0:0])

        interval_events = []

        for _, row in predicted_sequence.iterrows():
            start, end, class_id = int(row["start"]), int(row["end"]), int(row["class"])
            if start >= end:
                raise ValueError(
                    "Invalid predicted interval for "
                    f"sequence '{row['name']}', class {class_id}: start={start}, end={end}"
                )
            interval_events.append((start, "predicted_start", class_id, row.to_dict()))
            interval_events.append((end, "predicted_end", class_id, row.to_dict()))

        for _, row in true_sequence.iterrows():
            start, end, class_id = int(row["start"]), int(row["end"]), int(row["class"])
            if start >= end:
                raise ValueError(
                    "Invalid true interval for "
                    f"sequence '{row['name']}', class {class_id}: start={start}, end={end}"
                )
            interval_events.append((start, "true_start", class_id, row.to_dict()))
            interval_events.append((end, "true_end", class_id, row.to_dict()))

        interval_events.sort(key=lambda x: (x[0], EVENT_PRIORITY[x[1]]))

        # Atom intervals are assumed not to overlap within a sequence, so one active row per class is sufficient.
        active_predicted_classes = {}
        active_true_classes = {}
        previous_position = None

        for current_position, event_type, class_id, row_data in interval_events:
            if previous_position is not None and current_position > previous_position:
                interval_length = current_position - previous_position

                common_classes = set(active_predicted_classes.keys()).intersection(active_true_classes.keys())
                for atom_class in common_classes:
                    tp[atom_class] = tp.get(atom_class, 0) + interval_length
                    predicted_intervals.append(
                        _create_new_row(
                            row=active_predicted_classes[atom_class],
                            start=previous_position,
                            end=current_position,
                            status="TP"
                        )
                    )
                    true_intervals.append(
                        _create_new_row(
                            row=active_true_classes[atom_class],
                            start=previous_position,
                            end=current_position,
                            status="TP"
                        )
                    )

                for atom_class in set(active_predicted_classes.keys()) - common_classes:
                    fp[atom_class] = fp.get(atom_class, 0) + interval_length
                    predicted_intervals.append(
                        _create_new_row(
                            row=active_predicted_classes[atom_class],
                            start=previous_position,
                            end=current_position,
                            status="FP"
                        )
                    )

                for atom_class in set(active_true_classes.keys()) - common_classes:
                    fn[atom_class] = fn.get(atom_class, 0) + interval_length
                    true_intervals.append(
                        _create_new_row(
                            row=active_true_classes[atom_class],
                            start=previous_position,
                            end=current_position,
                            status="FN"
                        )
                    )

            if event_type == "predicted_start":
                active_predicted_classes[class_id] = row_data
            elif event_type == "predicted_end":
                active_predicted_classes.pop(class_id, None)
            elif event_type == "true_start":
                active_true_classes[class_id] = row_data
            elif event_type == "true_end":
                active_true_classes.pop(class_id, None)

            previous_position = current_position

    predicted_file = output_directory / "base_predicted_statuses.tsv"
    true_file = output_directory / "base_true_statuses.tsv"

    # -------- predicted intervals --------
    predicted_out = pd.DataFrame(predicted_intervals)
    _write_metrics_tsv(df=predicted_out, output_path=predicted_file)

    # -------- true intervals --------
    true_out = pd.DataFrame(true_intervals)
    _write_metrics_tsv(df=true_out, output_path=true_file)

    return tp, fp, fn
