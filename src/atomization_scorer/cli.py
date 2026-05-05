#!/usr/bin/env python3
"""
cli.py

Command-line interface (CLI) for Atomization Scorer.

Validates input genomes FASTA and GEESE atomization files.
Automatically creates an output directory if missing.
Computes overall atomization score using core scoring functions.

Functions
---------
validate_file       : Validate that a file exists and optionally check its extension.
validate_directory  : Validate that a directory exists; create if missing.
main                : CLI entry point, parses arguments, validates inputs, calls scoring.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
from __future__ import annotations
import argparse
import datetime
import json
import logging
import sys
from importlib.metadata import version as _package_version
from pathlib import Path

from atomization_scorer import compute_overall_score

# --------------------------------------------------------------------------------------
# Logging
# --------------------------------------------------------------------------------------
log = logging.getLogger(__name__)


# --------------------------------------------------------------------------------------
# Validation Functions
# --------------------------------------------------------------------------------------
def validate_file(path: Path, description: str, extension: str | tuple[str, ...]) -> None:
    """
    Validates that the input file exists and optionally has the required extension.

    Parameters
    ----------
    path : Path
        Path to the file to validate.
    description : str
        Name of the file for error messages.
    extension : str or tuple of str, optional, default=None
        Required file extension(s) (e.g., ".fa" or (".fa", ".fasta")), default is None.

    Raises
    ------
    SystemExit
        Exits program with sys.exit(1) if a file does not exist or the extension is incorrect.
    """
    if not path.is_file():
        log.error("%s file not found: %s", description, path)
        sys.exit(1)

    if extension:
        if isinstance(extension, str):
            extensions = (extension,)
        else:
            extensions = extension
        if path.suffix not in extensions:
            allowed = ", ".join(extensions)
            log.error("%s file must have one of the following extensions %s: %s", description, allowed, path)
            sys.exit(1)


def validate_directory(path: Path) -> None:
    """
    Validates that the output directory exists, creates it if missing.

    Parameters
    ----------
    path : Path
        Path to the output directory.

    Notes
    -----
    Prints a warning if the directory does not exist and creates it automatically.
    """
    if not path.is_dir():
        log.warning("Output directory '%s' does not exist. Creating directory...", path)
        path.mkdir(parents=True, exist_ok=True)
        log.info("Directory '%s' created.", path)


# --------------------------------------------------------------------------------------
# CLI Entry Point
# --------------------------------------------------------------------------------------
def main() -> None:
    """
    Parses command-line arguments, validates files/directories, and computes overall score.
    """
    parser = argparse.ArgumentParser(
        prog="atomization_scorer",
        description=(
            "Atomization Scorer - evaluate genome atomization quality.\n\n"
            "Computes an overall score from alignment-based and coverage-based metrics. "
            "The true (gold-standard) atomization is derived from representative-based "
            "minimap2 alignments. All pipeline parameters are tunable below."
        ),
        epilog=(
            "Examples:\n"
            "  # Basic usage with default settings\n"
            "  atomization_scorer genomes.fa predicted.geese ./results\n"
            "\n"
            "  # Base-level scoring with per-class breakdown\n"
            "  atomization_scorer genomes.fa predicted.geese ./results --level base --per-class\n"
            "\n"
            "  # Adjust weights (must sum to 1.0)\n"
            "  atomization_scorer genomes.fa predicted.geese ./results --alignment-weight 0.8 --coverage-weight 0.2\n"
            "\n"
            "  # Enable overlap diagnostics with a custom overlap threshold\n"
            "  atomization_scorer genomes.fa predicted.geese ./results --overlap-diagnostics --min-overlap-ratio 0.9\n"
            "\n"
            "  # Quiet output for use in pipelines\n"
            "  atomization_scorer genomes.fa predicted.geese ./results --quiet\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {_package_version('atomization_scorer')}",
    )

    verbosity_group = parser.add_mutually_exclusive_group()
    verbosity_group.add_argument(
        "--verbose",
        action="store_true",
        help="Enable DEBUG-level logging.",
    )
    verbosity_group.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress INFO messages; only WARNING and above are shown.",
    )

    # ------------------------------------------------------------------
    # Required positional arguments
    # ------------------------------------------------------------------
    parser.add_argument(
        "genomes_sequence",
        type=Path,
        metavar="GENOMES_FASTA",
        help="Input genomes FASTA file (.fa / .fasta).",
    )
    parser.add_argument(
        "geese_atomization",
        type=Path,
        metavar="ATOMIZATION_GEESE",
        help="Input GEESE atomization file (.geese).",
    )
    parser.add_argument(
        "output_directory",
        type=Path,
        metavar="OUTPUT_DIRECTORY",
        help=(
            "Output directory for results. Created if missing. "
            "Subdirectories: true_atomization/, metrics/, visualization/, "
            "overlap_diagnostics/ (if --overlap-diagnostics). "
            "Root files: scores.tsv, run_parameters.json."
        ),
    )

    # ------------------------------------------------------------------
    # Scoring options
    # ------------------------------------------------------------------
    scoring_group = parser.add_argument_group(
        "Scoring",
        "Control how alignment and coverage scores are computed and combined.",
    )
    scoring_group.add_argument(
        "--level",
        choices=["interval", "base"],
        default="interval",
        help=(
            "Alignment scoring granularity. "
            "'interval' matches whole predicted atoms to true atoms (default). "
            "'base' measures per-nucleotide overlap."
        ),
    )
    scoring_group.add_argument(
        "--per-class",
        action="store_true",
        help="Report alignment, coverage, and overall scores per atomization class instead of a single combined score.",
    )
    scoring_group.add_argument(
        "--min-overlap-ratio",
        dest="minimum_overlap_ratio",
        type=float,
        default=0.8,
        metavar="FLOAT",
        help=(
            "Minimum overlap ratio (intersection / union) for a predicted interval to count "
            "as a true positive in interval-level scoring. Default: 0.8."
        ),
    )
    scoring_group.add_argument(
        "--alignment-weight",
        type=float,
        default=0.7,
        metavar="FLOAT",
        help="Weight of the alignment score in the weighted geometric mean. Default: 0.7.",
    )
    scoring_group.add_argument(
        "--coverage-weight",
        type=float,
        default=0.3,
        metavar="FLOAT",
        help=(
            "Weight of the coverage score in the weighted geometric mean. "
            "Must satisfy alignment-weight + coverage-weight = 1.0. Default: 0.3."
        ),
    )

    # ------------------------------------------------------------------
    # True alignment pipeline options
    # ------------------------------------------------------------------
    pipeline_group = parser.add_argument_group(
        "True Alignment Pipeline",
        "Control representative selection and PAF filtering used to build the gold-standard atomization.",
    )
    pipeline_group.add_argument(
        "--representative-mode",
        choices=["mash", "first"],
        default="mash",
        dest="representative_mode",
        help=(
            "Strategy for picking one representative per predicted class. "
            "'mash' selects the atom with smallest total Mash distance to all others (default). "
            "'first' picks the first atom in the class."
        ),
    )
    pipeline_group.add_argument(
        "--min-similarity",
        type=float,
        default=0.95,
        metavar="FLOAT",
        dest="minimum_similarity",
        help="Minimum alignment similarity (0–1) to keep a PAF row after filtering. Default: 0.95.",
    )
    pipeline_group.add_argument(
        "--min-alignment-length",
        type=int,
        default=500,
        metavar="INT",
        dest="minimum_alignment_length",
        help="Minimum aligned-bases count to keep a PAF row after filtering. Default: 500.",
    )

    # ------------------------------------------------------------------
    # Minimap2 options
    # ------------------------------------------------------------------
    minimap2_group = parser.add_argument_group(
        "Minimap2",
        "Fine-tune the minimap2 alignment call.",
    )
    minimap2_group.add_argument(
        "--minimap2-preset",
        default="asm20",
        metavar="PRESET",
        dest="minimap2_preset",
        help="Minimap2 preset passed to -x (e.g. asm5, asm10, asm20, map-pb). Default: asm20.",
    )
    minimap2_group.add_argument(
        "--minimap2-secondary-ratio",
        type=float,
        default=0.1,
        metavar="FLOAT",
        dest="minimap2_secondary_ratio",
        help="Minimap2 secondary-to-primary score ratio passed to -p. Default: 0.1.",
    )
    minimap2_group.add_argument(
        "--minimap2-no-cigar",
        action="store_true",
        dest="minimap2_no_cigar",
        help="Omit CIGAR strings from PAF output (disables -c flag). Reduces file size but loses alignment detail.",
    )

    # ------------------------------------------------------------------
    # Overlap diagnostics options
    # ------------------------------------------------------------------
    diag_group = parser.add_argument_group(
        "Overlap Diagnostics",
        (
            "Generate diagnostic reports and dotplot inputs for overlapping alignments "
            "in the filtered PAF before overlap resolution. "
            "Outputs go to <output>/overlap_diagnostics/."
        ),
    )
    diag_group.add_argument(
        "--overlap-diagnostics",
        action="store_true",
        help="Enable overlap diagnostic reports and dotplot FASTA generation.",
    )
    diag_group.add_argument(
        "--overlap-report-min-len",
        type=int,
        default=0,
        metavar="INT",
        dest="minimum_report_overlap_length",
        help=(
            "Only include overlaps longer than this threshold (bp) in diagnostic reports. "
            "Default: 0 (report all)."
        ),
    )
    diag_group.add_argument(
        "--overlap-plot-min-len",
        type=int,
        default=0,
        metavar="INT",
        dest="minimum_plot_overlap_length",
        help=(
            "Only generate dotplot FASTA for overlaps longer than this threshold (bp). "
            "Default: 0 (plot all)."
        ),
    )
    diag_group.add_argument(
        "--overlap-include-reverse",
        action="store_true",
        dest="overlap_include_reverse",
        help="Also generate diagnostics from the partner's perspective (A→B and B→A).",
    )
    diag_group.add_argument(
        "--skip-dotter",
        action="store_true",
        dest="skip_dotter",
        help=(
            "Skip running Dotter (dot-plot visualization tool that requires Docker) "
            "after generating anchor FASTA inputs."
        ),
    )

    # Parse CLI arguments
    args = parser.parse_args()

    if args.verbose:
        log_level = logging.DEBUG
    elif args.quiet:
        log_level = logging.WARNING
    else:
        log_level = logging.INFO

    logging.basicConfig(
        level=log_level,
        format="%(asctime)s  %(levelname)-8s  %(message)s",
        datefmt="%H:%M:%S",
    )

    # ------------------------------------------------------------------
    # Validate inputs
    # ------------------------------------------------------------------
    validate_file(path=args.genomes_sequence, description="Genomes FASTA", extension=(".fa", ".fasta"))
    validate_file(path=args.geese_atomization, description="GEESE atomization", extension=".geese")
    validate_directory(path=args.output_directory)

    # ------------------------------------------------------------------
    # Info to user
    # ------------------------------------------------------------------
    log.info(
        "Processing files: genomes=%s atomization=%s output=%s",
        args.genomes_sequence,
        args.geese_atomization,
        args.output_directory,
    )

    # ------------------------------------------------------------------
    # Write run parameters for reproducibility
    # ------------------------------------------------------------------
    run_parameters = {
        "timestamp": datetime.datetime.now().isoformat(),
        "inputs": {
            "genomes_file": str(args.genomes_sequence.resolve()),
            "atomization_file": str(args.geese_atomization.resolve()),
        },
        "scoring": {
            "level": args.level,
            "per_class": args.per_class,
            "minimum_overlap_ratio": args.minimum_overlap_ratio,
            "alignment_weight": args.alignment_weight,
            "coverage_weight": args.coverage_weight,
        },
        "pipeline": {
            "representative_mode": args.representative_mode,
            "minimum_similarity": args.minimum_similarity,
            "minimum_alignment_length": args.minimum_alignment_length,
        },
        "minimap2": {
            "preset": args.minimap2_preset,
            "secondary_ratio": args.minimap2_secondary_ratio,
            "emit_cigar": not args.minimap2_no_cigar,
        },
        "diagnostics": {
            "run_overlap_diagnostics": args.overlap_diagnostics,
            "minimum_report_overlap_length": args.minimum_report_overlap_length,
            "minimum_plot_overlap_length": args.minimum_plot_overlap_length,
            "overlap_include_reverse": args.overlap_include_reverse,
            "run_dotter": not args.skip_dotter,
        },
    }
    run_parameters_file = args.output_directory / "run_parameters.json"
    with run_parameters_file.open("w") as file:
        json.dump(run_parameters, file, indent=2)
        file.write("\n")

    # ------------------------------------------------------------------
    # Compute overall score
    # ------------------------------------------------------------------
    log.info("=" * 60)
    log.info("Computing overall score...")
    result = compute_overall_score(
        genomes_file=args.genomes_sequence,
        atomization_file=args.geese_atomization,
        output_directory=args.output_directory,
        level=args.level,
        per_class=args.per_class,
        minimum_overlap_ratio=args.minimum_overlap_ratio,
        alignment_weight=args.alignment_weight,
        coverage_weight=args.coverage_weight,
        representative_mode=args.representative_mode,
        minimum_similarity=args.minimum_similarity,
        minimum_alignment_length=args.minimum_alignment_length,
        minimap2_preset=args.minimap2_preset,
        minimap2_secondary_ratio=args.minimap2_secondary_ratio,
        minimap2_emit_cigar=not args.minimap2_no_cigar,
        run_overlap_diagnostics=args.overlap_diagnostics,
        minimum_report_overlap_length=args.minimum_report_overlap_length,
        minimum_plot_overlap_length=args.minimum_plot_overlap_length,
        overlap_include_reverse=args.overlap_include_reverse,
        run_dotter=not args.skip_dotter,
    )

    log.info("=" * 60)
    log.info("Overall score result: %s", result)


# --------------------------------------------------------------------------------------
# Execute CLI if script is run directly
# --------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
