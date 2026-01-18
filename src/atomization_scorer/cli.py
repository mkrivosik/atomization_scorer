#!/usr/bin/env python3
"""
cli.py

Command-line interface (CLI) for Atomization Scorer.

Validates input genomes FASTA and GEESE atomization files.
Automatically creates output directory if missing.
Computes overall atomization score using core scoring functions.

Functions
---------
validate_file       : Validate that a file exists and optionally check its extension.
validate_directory  : Validate that a directory exists; create if missing.
main                : CLI entry point, parses arguments, validates inputs, calls scoring.
"""

# ---------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------
import argparse
import sys
from typing import Tuple
from pathlib import Path
from atomization_scorer import compute_overall_score

# ---------------------------------------------------------------------
# Validation Functions
# ---------------------------------------------------------------------

def validate_file(path: Path, description: str, extension: str | Tuple[str, ...]) -> None:
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
        Exits program with sys.exit(1) if file does not exist or extension is incorrect.
    """
    if not path.is_file():
        print(f"Error: {description} file not found: {path}", file=sys.stderr)
        sys.exit(1)

    if extension:
        if isinstance(extension, str):
            extensions = (extension,)
        else:
            extensions = extension
        if path.suffix not in extensions:
            allowed = ", ".join(extensions)
            print(f"Error: {description} file must have one of the following extensions {allowed}: {path}", file=sys.stderr)
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
    Prints a warning if directory does not exist and creates it automatically.
    """
    if not path.is_dir():
        print(f"Warning: Output directory '{path}' does not exist. Creating directory...", file=sys.stderr)
        path.mkdir(parents=True, exist_ok=True)
        print(f"Directory '{path}' created.", file=sys.stderr)

# ---------------------------------------------------------------------
# CLI Entry Point
# ---------------------------------------------------------------------

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

    # ---------------------------
    # Validate inputs
    # ---------------------------
    validate_file(path=args.genomes_sequence, description="Genomes FASTA", extension=(".fa", ".fasta"))
    validate_file(path=args.geese_atomization, description="GEESE atomization", extension=".geese")
    validate_directory(path=args.output_directory)

    # ---------------------------
    # Info to user
    # ---------------------------
    print("Processing files:")
    print(f"  Genomes FASTA:        {args.genomes_sequence}")
    print(f"  GEESE atomization:    {args.geese_atomization}")
    print(f"  Output directory:     {args.output_directory}")

    # ---------------------------
    # Compute overall score
    # ---------------------------
    print("\nComputing overall score...")
    result = compute_overall_score(
        genomes_file=args.genomes_sequence,
        atomization_file=args.geese_atomization,
        output_directory=args.output_directory
    )

    print("Overall score result: ", result)

# ---------------------------------------------------------------------
# Execute CLI if script is run directly
# ---------------------------------------------------------------------

if __name__ == "__main__":
    main()
