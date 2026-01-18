"""
pipeline

This subpackage provides the alignment pipelines for genome atomization scoring.

Modules
-------
true_pipeline               : Computes gold-standard (true) genome atomization alignment.
minimap2_aligner            : Computes minimap2 alignment.
representatives_selector    : Extracts representative sequences.
"""

# ---------------------------------------------------------------------
# True Pipeline
# ---------------------------------------------------------------------
from .true_pipeline import compute_true_alignment

# ---------------------------------------------------------------------
# Minimap2 Aligner
# ---------------------------------------------------------------------
from .minimap2_aligner import align_with_minimap2

# ---------------------------------------------------------------------
# Representatives Selector
# ---------------------------------------------------------------------
from .representatives_selector import extract_representatives

# ---------------------------------------------------------------------
# Package API
# ---------------------------------------------------------------------
__all__ = [
    # True pipeline
    'compute_true_alignment',

    # Minimap2 aligner
    'align_with_minimap2',

    # Representatives selector
    'extract_representatives',
]
