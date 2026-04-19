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
import logging
import sys
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

    Command-line arguments
    ----------------------
    genomes_sequence    : Input genomes FASTA file (.fa/.fasta).
    geese_atomization   : Input GEESE atomization file (.geese).
    output_directory    : Directory where results will be stored.
    """
    parser = argparse.ArgumentParser(
        prog="atomization_scorer",
        description="Atomization Scorer - tool for evaluating genomes atomization."
    )

    # Define CLI arguments
    parser.add_argument("genomes_sequence", type=Path, help="Input genomes FASTA file (.fa/.fasta).")
    parser.add_argument("geese_atomization", type=Path, help="Input GEESE atomization file (.geese).")
    parser.add_argument("output_directory", type=Path, help="Output directory for results.")

    # Parse CLI arguments
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-8s  %(message)s",
        datefmt="%H:%M:%S",
    )

    # ---------------------------
    # Validate inputs
    # ---------------------------
    validate_file(path=args.genomes_sequence, description="Genomes FASTA", extension=(".fa", ".fasta"))
    validate_file(path=args.geese_atomization, description="GEESE atomization", extension=".geese")
    validate_directory(path=args.output_directory)

    # ---------------------------
    # Info to user
    # ---------------------------
    log.info(
        "Processing files: genomes=%s atomization=%s output=%s",
        args.genomes_sequence,
        args.geese_atomization,
        args.output_directory,
    )

    # ---------------------------
    # Compute overall score
    # ---------------------------
    log.info("=" * 60)
    log.info("Computing overall score...")
    result = compute_overall_score(
        genomes_file=args.genomes_sequence,
        atomization_file=args.geese_atomization,
        output_directory=args.output_directory
    )

    log.info("=" * 60)
    log.info("Overall score result: %s", result)


# --------------------------------------------------------------------------------------
# Execute CLI if script is run directly
# --------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
