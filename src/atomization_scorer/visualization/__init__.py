"""
visualization

This package provides functions for visualizing genome atomization results.

Modules
-------
genome_visualization : Function for plotting genome atomization intervals.
"""

# ---------------------------------------------------------------------
# Genome Visualization
# ---------------------------------------------------------------------
from .atomization_visualization import plot_genome_atomization

# ---------------------------------------------------------------------
# Package API
# ---------------------------------------------------------------------
__all__ = [
    'plot_genome_atomization',
]
