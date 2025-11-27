"""
Guide extraction and perturbation assignment for Crop-Seq.

This module provides functions to identify CRISPR guide RNAs in single cells
and assign perturbations based on guide identity.
"""

from perturbio.guides.extraction import extract, assign_guides_to_cells, detect_multiplets

__all__ = [
    "extract",
    "assign_guides_to_cells",
    "detect_multiplets",
]
