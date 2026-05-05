"""
data_processing

This package provides functions for reading and processing input data files
used in genome atomization scoring.

Modules
-------
fasta_reader     : Module for reading FASTA genome sequences.
fasta_writer     : Module for writing FASTA sequence files.
geese_reader     : Module for reading GEESE atomization files.
geese_writer     : Module for writing GEESE atomization files.
paf_processing   : Module for filtering Minimap2 PAF files.
paf_segmentation : Module for resolving PAF overlaps and validating true atomization.
paf_to_geese     : Module for converting PAF files to GEESE format.
utils            : Module providing shared utility functions.
"""

# --------------------------------------------------------------------------------------
# FASTA Reader
# --------------------------------------------------------------------------------------
from .fasta_reader import read_fasta

# --------------------------------------------------------------------------------------
# FASTA Writer
# --------------------------------------------------------------------------------------
from .fasta_writer import write_fasta

# --------------------------------------------------------------------------------------
# GEESE Reader
# --------------------------------------------------------------------------------------
from .geese_reader import read_geese

# --------------------------------------------------------------------------------------
# GEESE Writer
# --------------------------------------------------------------------------------------
from .geese_writer import write_geese

# --------------------------------------------------------------------------------------
# PAF Processing
# --------------------------------------------------------------------------------------
from .paf_processing import filter_paf

# --------------------------------------------------------------------------------------
# PAF To GEESE
# --------------------------------------------------------------------------------------
from .paf_to_geese import paf_to_geese

# --------------------------------------------------------------------------------------
# PAF Segmentation
# --------------------------------------------------------------------------------------
from .paf_segmentation import (
    resolve_paf_overlaps,
    validate_non_overlapping_geese,
)

# --------------------------------------------------------------------------------------
# Utils
# --------------------------------------------------------------------------------------
from .utils import (
    check_required_columns,
)

# --------------------------------------------------------------------------------------
# Package API
# --------------------------------------------------------------------------------------
__all__ = [
    # FASTA reader
    'read_fasta',

    # FASTA writer
    'write_fasta',

    # GEESE reader
    'read_geese',

    # GEESE writer
    'write_geese',

    # PAF processing
    'filter_paf',

    # PAF segmentation
    'resolve_paf_overlaps',
    'validate_non_overlapping_geese',

    # PAF processing
    'paf_to_geese',

    # Utility functions
    'check_required_columns',
]
