"""
scoring_system

This package provides functions for score calculations.

Modules
-------
alignment_score        : Module for computing alignment scores.
base_metrics           : Module for computing base-level metrics.
coverage_score         : Module for computing genome coverage score.
interval_metrics       : Module for computing interval-level metrics.
overall_score          : Module for computing overall atomization score.
artificial_degradation : Module for creating artificially degraded atomization.
"""

# --------------------------------------------------------------------------------------
# Alignment Score
# --------------------------------------------------------------------------------------
from .alignment_score import compute_alignment_score

# --------------------------------------------------------------------------------------
# Base Metrics
# --------------------------------------------------------------------------------------
from .base_metrics import compute_base_level_metrics

# --------------------------------------------------------------------------------------
# Coverage Score
# --------------------------------------------------------------------------------------
from .coverage_score import compute_coverage_score

# --------------------------------------------------------------------------------------
# Artificial Degradation
# --------------------------------------------------------------------------------------
from .artificial_degradation import degrade_atomization

# --------------------------------------------------------------------------------------
# Interval Metrics
# --------------------------------------------------------------------------------------
from .interval_metrics import compute_interval_level_metrics

# --------------------------------------------------------------------------------------
# Overall Score
# --------------------------------------------------------------------------------------
from .overall_score import compute_overall_score

# --------------------------------------------------------------------------------------
# Package API
# --------------------------------------------------------------------------------------
__all__ = [
    # Base metrics
    'compute_base_level_metrics',

    # Interval metrics
    'compute_interval_level_metrics',

    # Alignment score
    'compute_alignment_score',

    # Coverage score
    'compute_coverage_score',

    # Artificial degradation
    'degrade_atomization',

    # Overall score
    'compute_overall_score',
]
