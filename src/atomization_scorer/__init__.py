"""
atomization_scorer

This package provides the main scoring functions, alignment pipelines,
and data processing utilities for evaluating genome atomization.

Modules
-------
data_processing : Subpackage providing data readers and writers for PSL and GEESE files.
pipeline        : Subpackage providing predicted and true alignment pipelines.
scoring_system  : Subpackage providing scoring functions.
visualization   : Subpackage providing a genome atomization plotting function.
"""

# --------------------------------------------------------------------------------------
# Data Processing
# --------------------------------------------------------------------------------------
from .data_processing import (
    filter_paf,
    paf_to_geese,
    resolve_paf_overlaps,
    read_fasta,
    read_geese,
    validate_non_overlapping_geese,
    write_fasta,
)

# --------------------------------------------------------------------------------------
# Pipeline
# --------------------------------------------------------------------------------------
from .pipeline import (
    align_with_minimap2,
    compute_true_alignment,
    extract_representatives,
)

# --------------------------------------------------------------------------------------
# Scoring System
# --------------------------------------------------------------------------------------
from .scoring_system import (
    compute_alignment_score,
    compute_base_level_metrics,
    compute_coverage_score,
    compute_interval_level_metrics,
    compute_overall_score,
)

# --------------------------------------------------------------------------------------
# Visualization
# --------------------------------------------------------------------------------------
from .visualization import (
    plot_atomization,
)

# --------------------------------------------------------------------------------------
# Package API
# --------------------------------------------------------------------------------------
__all__ = [
    # Data processing
    'read_fasta',
    'write_fasta',
    'read_geese',
    'filter_paf',
    'resolve_paf_overlaps',
    'paf_to_geese',
    'validate_non_overlapping_geese',

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

    # Visualization
    'plot_atomization',
]
