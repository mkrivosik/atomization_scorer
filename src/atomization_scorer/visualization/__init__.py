"""
visualization

This package provides an interactive HTML visualization for genome atomization results.

Modules
-------
atomization_visualization : Module for rendering genome atomization plots.
"""

# --------------------------------------------------------------------------------------
# Visualization API
# --------------------------------------------------------------------------------------
from .atomization_visualization import plot_atomization

# --------------------------------------------------------------------------------------
# Package API
# --------------------------------------------------------------------------------------
__all__ = [
    'plot_atomization',
]
