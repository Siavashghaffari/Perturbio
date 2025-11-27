"""
Analysis tools for Crop-Seq data.

This module provides statistical analysis functions for identifying
genes affected by CRISPR perturbations.
"""

from perturbio.analysis.differential import differential_expression, rank_genes_per_perturbation

__all__ = [
    "differential_expression",
    "rank_genes_per_perturbation",
]
