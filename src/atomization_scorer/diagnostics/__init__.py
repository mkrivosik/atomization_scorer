"""
diagnostics

This subpackage provides diagnostic utilities for inspecting suspicious
true-atomization overlaps before overlap resolution.

Modules
-------
overlap_diagnostics : Module for generating overlap reports and dotplot FASTA inputs.
dotplot_inputs      : Module for building per-anchor FASTA files for dotplot rendering.
dotter_runner       : Module for running Dotter on per-anchor FASTA inputs.
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
