"""
pipeline

This subpackage provides the alignment pipelines for genome atomization scoring.

Modules
-------
true_pipeline             : Module for computing gold-standard true atomization.
minimap2_aligner          : Module for aligning genomes with minimap2.
representatives_selector  : Module for extracting representative sequences.
"""

# --------------------------------------------------------------------------------------
# Minimap2 Aligner
# --------------------------------------------------------------------------------------
from .minimap2_aligner import align_with_minimap2

# --------------------------------------------------------------------------------------
# Representatives Selector
# --------------------------------------------------------------------------------------
from .representatives_selector import extract_representatives

# --------------------------------------------------------------------------------------
# True Pipeline
# --------------------------------------------------------------------------------------
from .true_pipeline import compute_true_alignment

# --------------------------------------------------------------------------------------
# Package API
# --------------------------------------------------------------------------------------
__all__ = [
    # True pipeline
    'compute_true_alignment',

    # Minimap2 aligner
    'align_with_minimap2',

    # Representatives selector
    'extract_representatives',
]
