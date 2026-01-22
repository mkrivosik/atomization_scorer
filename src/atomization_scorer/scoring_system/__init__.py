"""
scoring_system

This package provides functions for score calculations.

Modules
-------
compute_base_level_metrics      : Function for computing base level metrics.
compute_interval_level_metrics  : Function for computing interval level metrics.
compute_alignment_score         : Function for computing alignment scores.
compute_coverage_score          : Function for computing genomes coverage score.
overall_score                   : Function for computing overall atomization score.
"""

# ---------------------------------------------------------------------
# Alignment Score
# ---------------------------------------------------------------------
from .alignment_score import compute_alignment_score

# ---------------------------------------------------------------------
# Base Metrics
# ---------------------------------------------------------------------
from .base_metrics import compute_base_level_metrics

# ---------------------------------------------------------------------
# Coverage Score
# ---------------------------------------------------------------------
from .coverage_score import compute_coverage_score

# ---------------------------------------------------------------------
# Interval Metrics
# ---------------------------------------------------------------------
from .interval_metrics import compute_interval_level_metrics

# ---------------------------------------------------------------------
# Overall Score
# ---------------------------------------------------------------------
from .overall_score import compute_overall_score

# ---------------------------------------------------------------------
# Package API
# ---------------------------------------------------------------------
__all__ = [
    # Base metrics
    'compute_base_level_metrics',

    # Interval metrics
    'compute_interval_level_metrics',

    # Alignment score
    'compute_alignment_score',

    # Coverage score
    'compute_coverage_score',

    # Overall score
    'compute_overall_score',
]
