"""
Log-scaled atom length distribution plot with good vs global singleton overlay.
"""

# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------
from __future__ import annotations

import argparse
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd

from atomization_scorer.data_processing import read_geese

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
HERE = Path(__file__).parent
FIXTURES = HERE.parent.parent.parent / "tests" / "fixtures"

DEFAULT_GEESE = FIXTURES / "big.geese"
DEFAULT_OUTPUT = HERE.parent / "outputs" / "atom_length_distribution.png"
DEFAULT_N_BINS = 60
TAIL_START_BP = 20_000

COLOR_ALL = "#4C72B0"
COLOR_MULTITON = "#4C72B0"
COLOR_SINGLETON = "#C44E52"
ALPHA_OVERLAY = 0.85
BAR_RWIDTH = 0.88

_RC = {
    "font.family": "sans-serif",
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.8,
    "axes.labelsize": 11,
    "axes.titlesize": 12,
    "axes.titleweight": "bold",
    "xtick.major.size": 4,
    "ytick.major.size": 4,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "legend.frameon": True,
    "legend.framealpha": 0.92,
    "legend.edgecolor": "#cccccc",
    "legend.fontsize": 9,
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _compute_lengths(atoms_df: pd.DataFrame) -> pd.DataFrame:
    """Add 'length' column computed as end - start."""
    df = atoms_df.copy()
    df["length"] = df["end"].astype(int) - df["start"].astype(int)
    return df


def _classify_atoms(atoms_df: pd.DataFrame) -> pd.DataFrame:
    """Add 'category' column: 'singleton' if class has exactly 1 atom globally, else 'multiton'."""
    class_counts = atoms_df.groupby("class").size()
    singleton_classes = set(class_counts[class_counts == 1].index)
    atoms_df = atoms_df.copy()
    atoms_df["category"] = atoms_df["class"].apply(
        lambda c: "singleton" if c in singleton_classes else "multiton"
    )
    return atoms_df


def _log_bins(lengths: np.ndarray, n_bins: int) -> np.ndarray:
    """Return n_bins log-spaced bin edges covering the range of lengths."""
    lower = np.log10(lengths.min())
    upper = np.log10(lengths.max())
    if lower == upper:
        return np.array([lengths.min(), lengths.max() + 1])
    return np.logspace(lower, upper, n_bins + 1)


def _bp_formatter(x: float, _pos: int) -> str:
    """Format a base-pair count as human-readable bp / kb / Mb."""
    if x >= 1_000_000:
        return f"{x / 1_000_000:.3g} Mb"
    if x >= 1_000:
        return f"{x / 1_000:.3g} kb"
    return f"{x:.3g} bp"


def _apply_log_x(ax: plt.Axes) -> None:
    """Set log scale on x-axis with human-readable bp/kb/Mb tick labels."""
    ax.set_xscale("log")
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(_bp_formatter))
    ax.tick_params(axis="x", labelrotation=30, labelsize=9)


def _apply_count_y(ax: plt.Axes) -> None:
    """Set y-axis label, comma-formatted tick labels, and a light horizontal grid."""
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
    ax.set_ylabel("Count")
    ax.grid(axis="y", alpha=0.25, linewidth=0.5, zorder=0)


def _panel_label(ax: plt.Axes, label: str) -> None:
    """Place a bold panel label (A, B, C, D) in the upper-left corner of the axes."""
    ax.text(-0.13, 1.04, label, transform=ax.transAxes,
            fontsize=14, fontweight="bold", va="bottom", ha="left")


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def _plot_full(ax: plt.Axes, lengths: np.ndarray, bins: np.ndarray) -> None:
    """Panel A: full atom length distribution over all atoms on a log x-axis."""
    ax.hist(lengths, bins=bins, color=COLOR_ALL, edgecolor="white",
            linewidth=0.25, rwidth=BAR_RWIDTH)
    _apply_log_x(ax)
    _apply_count_y(ax)
    ax.set_xlabel("Atom length")
    ax.set_title("Atom Length Distribution")

    stats = (
        f"$n$ = {len(lengths):,}\n"
        f"min = {_bp_formatter(lengths.min(), 0)}\n"
        f"median = {_bp_formatter(float(np.median(lengths)), 0)}\n"
        f"max = {_bp_formatter(lengths.max(), 0)}"
    )
    ax.text(0.97, 0.97, stats, transform=ax.transAxes,
            ha="right", va="top", fontsize=8.5, linespacing=1.55,
            bbox=dict(boxstyle="round,pad=0.45", facecolor="white",
                      alpha=0.92, edgecolor="#cccccc", linewidth=0.6))
    _panel_label(ax, "A")


def _plot_overlay(ax: plt.Axes, df: pd.DataFrame, bins: np.ndarray) -> None:
    """Panel B: overlaid multiton and singleton count histograms on a log x-axis."""
    multiton = df[df["category"] == "multiton"]["length"].to_numpy()
    singleton = df[df["category"] == "singleton"]["length"].to_numpy()
    total = len(df)

    ax.hist(multiton, bins=bins, alpha=ALPHA_OVERLAY, color=COLOR_MULTITON,
            edgecolor="white", linewidth=0.25, rwidth=BAR_RWIDTH,
            label=f"Multiton  ($n$={len(multiton):,},  {100 * len(multiton) / total:.1f}%)")
    ax.hist(singleton, bins=bins, alpha=ALPHA_OVERLAY, color=COLOR_SINGLETON,
            edgecolor="white", linewidth=0.25, rwidth=BAR_RWIDTH,
            label=f"Singleton  ($n$={len(singleton):,},  {100 * len(singleton) / total:.1f}%)")

    _apply_log_x(ax)
    _apply_count_y(ax)
    ax.set_xlabel("Atom length")
    ax.set_title("Multiton vs Singleton")
    ax.legend(loc="upper right")
    _panel_label(ax, "B")


def _plot_tail_zoom(ax: plt.Axes, df: pd.DataFrame) -> None:
    """Panel C: zoomed log x-axis view of atoms at or above TAIL_START_BP."""
    tail = df[df["length"] >= TAIL_START_BP].copy()
    multiton = tail[tail["category"] == "multiton"]["length"].to_numpy()
    singleton = tail[tail["category"] == "singleton"]["length"].to_numpy()
    bins = _log_bins(tail["length"].to_numpy(), n_bins=30)

    ax.hist(multiton, bins=bins, alpha=ALPHA_OVERLAY, color=COLOR_MULTITON,
            edgecolor="white", linewidth=0.25, rwidth=BAR_RWIDTH,
            label=f"Multiton  ($n$={len(multiton):,})")
    ax.hist(singleton, bins=bins, alpha=ALPHA_OVERLAY, color=COLOR_SINGLETON,
            edgecolor="white", linewidth=0.25, rwidth=BAR_RWIDTH,
            label=f"Singleton  ($n$={len(singleton):,})")

    ax.set_xscale("log")
    ax.set_xticks([20_000, 50_000, 100_000])
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(_bp_formatter))
    ax.tick_params(axis="x", labelrotation=30, labelsize=9)
    _apply_count_y(ax)
    ax.set_xlabel("Atom length")
    ax.set_title(f"Tail Region  (≥ {_bp_formatter(TAIL_START_BP, 0)})")
    ax.legend(loc="upper right")
    _panel_label(ax, "C")


def _plot_singleton_fraction(ax: plt.Axes, df: pd.DataFrame, bins: np.ndarray) -> None:
    """Panel D: singleton fraction per length bin with overall average as a reference line."""
    multiton_counts, _ = np.histogram(
        df[df["category"] == "multiton"]["length"].to_numpy(), bins=bins
    )
    singleton_counts, _ = np.histogram(
        df[df["category"] == "singleton"]["length"].to_numpy(), bins=bins
    )
    total = multiton_counts + singleton_counts
    safe_total = np.where(total > 0, total, 1)
    fraction = np.where(total > 0, singleton_counts / safe_total * 100, np.nan)

    centers = np.sqrt(bins[:-1] * bins[1:])
    widths = np.diff(bins) * BAR_RWIDTH
    ax.bar(centers, fraction, width=widths,
           color=COLOR_SINGLETON, alpha=0.8, edgecolor="white", linewidth=0.25)

    overall = (df["category"] == "singleton").mean() * 100
    ax.axhline(overall, color="#555555", linewidth=1.0, linestyle="--",
               label=f"Overall  ({overall:.1f}%)")

    _apply_log_x(ax)
    ax.set_xlabel("Atom length")
    ax.set_ylabel("Singleton fraction (%)")
    ax.set_title("Singleton Fraction by Length")
    ax.set_ylim(bottom=0)
    ax.grid(axis="y", alpha=0.25, linewidth=0.5, zorder=0)
    ax.legend(loc="upper right")
    _panel_label(ax, "D")


def plot_distribution(atoms_df: pd.DataFrame, output: Path, n_bins: int = DEFAULT_N_BINS) -> None:
    """Render the 2×2 panel figure and save it as a PNG."""
    df = _compute_lengths(atoms_df)
    df = _classify_atoms(df)
    lengths = df["length"].to_numpy()
    bins = _log_bins(lengths, n_bins)

    with plt.rc_context(_RC):
        fig, axes = plt.subplots(2, 2, figsize=(13, 8))
        fig.suptitle("Atom Length Distributions", fontsize=13, fontweight="bold", y=0.99)
        fig.subplots_adjust(wspace=0.38, hspace=0.55, left=0.08, right=0.97, top=0.91, bottom=0.12)

        _plot_full(axes[0, 0], lengths, bins)
        _plot_overlay(axes[0, 1], df, bins)
        _plot_tail_zoom(axes[1, 0], df)
        _plot_singleton_fraction(axes[1, 1], df, bins)

        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=300, bbox_inches="tight")
        plt.close(fig)
    log.info("Plot saved to %s", output)


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
def run(geese: Path, output: Path, n_bins: int = DEFAULT_N_BINS) -> None:
    """Load atoms from a GEESE file and produce the distribution figure."""
    log.info("Reading atoms from %s", geese)
    atoms_df = read_geese(geese_file=geese)
    log.info("Loaded %d atoms", len(atoms_df))
    plot_distribution(atoms_df, output=output, n_bins=n_bins)


# ---------------------------------------------------------------------------
# Command Line Interface
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Plot atom length distribution (log scale).")
    parser.add_argument("--geese", type=Path, default=DEFAULT_GEESE, help="Atoms GEESE file.")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="Output PNG path.")
    parser.add_argument("--n-bins", type=int, default=DEFAULT_N_BINS, help="Number of histogram bins.")
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Script Entry Point
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-8s  %(message)s",
        datefmt="%H:%M:%S",
    )
    args = parse_args()
    run(geese=args.geese, output=args.output, n_bins=args.n_bins)
