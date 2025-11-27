"""
Perturbio: Comprehensive Crop-Seq Analysis in Python
=====================================================

Perturbio provides tools for end-to-end analysis of Crop-Seq experiments
(CRISPR pooled screens + single-cell RNA sequencing).

Main Features
-------------
- Guide extraction from single-cell data
- Differential expression analysis
- Publication-quality visualizations
- Scanpy-compatible workflows

Quick Start
-----------
>>> from perturbio import CropSeqAnalyzer
>>> analyzer = CropSeqAnalyzer("cropseq_data.h5ad")
>>> results = analyzer.run()
>>> results.plot_volcano("BRCA1_guide1")

Scanpy Integration
------------------
>>> import scanpy as sc
>>> import perturbio as pt
>>> adata = sc.read_h5ad("data.h5ad")
>>> pt.guides.extract(adata, guide_file="guides.csv")
>>> pt.tl.differential_expression(adata, groupby='perturbation')
"""

__version__ = "0.1.0"
__author__ = "Perturbio Developers"

# High-level API
from perturbio.core import CropSeqAnalyzer

# Submodules
from perturbio import io
from perturbio import guides
from perturbio import analysis as tl  # 'tl' for 'tools' (scanpy convention)
from perturbio import plotting as pl  # 'pl' for 'plotting' (scanpy convention)
from perturbio import results

# Public API
__all__ = [
    # Version
    "__version__",
    # High-level API
    "CropSeqAnalyzer",
    # Submodules
    "io",
    "guides",
    "tl",  # analysis tools
    "pl",  # plotting
    "results",
]
