"""
diagnostics

This subpackage provides diagnostic utilities for inspecting suspicious
true-atomization overlaps before overlap resolution.

Modules
-------
overlap_diagnostics : Generate overlap reports and dotplot FASTA inputs from filtered PAF files.
dotplot_inputs      : Build per-anchor FASTA files for downstream dotplot rendering.
dotter_runner       : Run Dotter on per-anchor FASTA inputs.
"""

# --------------------------------------------------------------------------------------
# Overlap Diagnostics
# --------------------------------------------------------------------------------------
from .overlap_diagnostics import diagnose_paf_overlaps

# --------------------------------------------------------------------------------------
# Dotter Runner
# --------------------------------------------------------------------------------------
from .dotter_runner import (
    run_dotter_for_anchor,
    run_dotter_for_anchors,
)

# --------------------------------------------------------------------------------------
# Package API
# --------------------------------------------------------------------------------------
__all__ = [
    'diagnose_paf_overlaps',
    'run_dotter_for_anchor',
    'run_dotter_for_anchors',
]
