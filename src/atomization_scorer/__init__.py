"""
atomization_scorer

This package provides the main scoring functions, alignment pipelines,
and data processing utilities for evaluating genome atomization.

Modules
-------
data_processing : Subpackage providing data readers and writers for PSL and GEESE files.
pipeline        : Subpackage providing predicted and true alignment pipelines.
scoring_system  : Subpackage providing scoring functions.
"""

# ---------------------------------------------------------------------
# Data Processing
# ---------------------------------------------------------------------
from .data_processing import (
    read_fasta,
    write_fasta,
    read_geese,
    filter_paf,
    paf_to_geese,
)

# ---------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------
from .pipeline import (
    compute_true_alignment,
    align_with_minimap2,
    extract_representatives,
)

# ---------------------------------------------------------------------
# Scoring System
# ---------------------------------------------------------------------
from .scoring_system import (
    compute_base_level_metrics,
    compute_interval_level_metrics,
    compute_alignment_score,
    compute_coverage_score,
    compute_overall_score,
)

# ---------------------------------------------------------------------
# Package API
# ---------------------------------------------------------------------
__all__ = [
    # Data processing
    'read_fasta',
    'write_fasta',
    'read_geese',
    'filter_paf',
    'paf_to_geese',

    # Pipeline
    'compute_true_alignment',
    'align_with_minimap2',
    'extract_representatives',

    # Scoring system
    'compute_base_level_metrics',
    'compute_interval_level_metrics',
    'compute_alignment_score',
    'compute_coverage_score',
    'compute_overall_score',
]
