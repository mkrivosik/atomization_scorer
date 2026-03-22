"""
visualization

This package provides a function for visualizing genome atomization results.

Modules
-------
atomization_visualization   : Function for plotting genome atomization.
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
