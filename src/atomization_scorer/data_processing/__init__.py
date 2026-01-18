"""
data_processing

This package provides functions for reading and processing input data files
used in genome atomization scoring.

Modules
-------
fasta_reader    : Function for reading FASTA genome sequences.
fasta_writer    : Function for writing FASTA sequence files.
geese_reader    : Function for reading GEESE atomization TSV files.
paf_processing  : Function for processing Minimap2 PAF files.
paf_to_geese    : Function for converting PAF files to GEESE files.
utils           : Shared utility functions (column checks, renaming, etc.)
"""

# ---------------------------------------------------------------------
# FASTA Reader
# ---------------------------------------------------------------------
from .fasta_reader import read_fasta

# ---------------------------------------------------------------------
# FASTA Writer
# ---------------------------------------------------------------------
from .fasta_writer import write_fasta

# ---------------------------------------------------------------------
# GEESE Reader
# ---------------------------------------------------------------------
from .geese_reader import read_geese

# ---------------------------------------------------------------------
# PAF Processing
# ---------------------------------------------------------------------
from .paf_processing import filter_paf

# ---------------------------------------------------------------------
# PAF To GEESE
# ---------------------------------------------------------------------
from .paf_to_geese import paf_to_geese

# ---------------------------------------------------------------------
# Utils
# ---------------------------------------------------------------------
from .utils import (
    check_required_columns,
    rename_column,
)

# ---------------------------------------------------------------------
# Package API
# ---------------------------------------------------------------------
__all__ = [
    # FASTA reader
    'read_fasta',

    # FASTA writer
    'write_fasta',

    # GEESE reader
    'read_geese',

    # PAF processing
    'filter_paf',

    # PAF processing
    'paf_to_geese',

    # Utility functions
    'check_required_columns',
    'rename_column',
]
