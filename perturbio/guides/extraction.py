"""
Guide extraction algorithm for Crop-Seq data.

This module implements the core algorithm for identifying which CRISPR guide
perturbed each cell based on guide RNA expression.
"""

from typing import Union, Optional
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import anndata as ad
from tqdm import tqdm

from perturbio.io.readers import load_guide_library


def extract(
    adata: ad.AnnData,
    guide_file: Optional[Union[str, Path]] = None,
    guide_library: Optional[pd.DataFrame] = None,
    min_umis: int = 3,
    multiplet_threshold: float = 0.2,
    key_added: str = "perturbation",
    copy: bool = False,
) -> Optional[ad.AnnData]:
    """
    Extract guide barcodes from cells and assign perturbations.

    This function identifies which CRISPR guide perturbed each cell by looking
    for guide RNA expression. It uses a winner-takes-all approach: the guide
    with the highest UMI count in each cell is assigned to that cell.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with cells × genes
    guide_file : str or Path, optional
        Path to guide library CSV file. Either this or guide_library must be provided.
    guide_library : DataFrame, optional
        Guide library DataFrame with columns: guide_id, target_gene
    min_umis : int, default 3
        Minimum UMI count required to confidently assign a guide
    multiplet_threshold : float, default 0.2
        Threshold for detecting multiplets. If a cell has 2+ guides with UMI counts
        within this fraction of the max count, it's flagged as a multiplet.
    key_added : str, default "perturbation"
        Key name for storing perturbation annotations in adata.obs
    copy : bool, default False
        Whether to return a copy of adata or modify in place

    Returns
    -------
    AnnData or None
        If copy=True, returns modified copy of adata.
        If copy=False, modifies adata in place and returns None.

    Notes
    -----
    Adds the following to adata.obs:
    - `{key_added}`: Target gene for the assigned guide
    - `guide_identity`: Guide ID assigned to each cell
    - `guide_umi_count`: UMI count for the assigned guide
    - `guide_confidence`: Assignment confidence (high, low, multiplet, unassigned)

    Adds to adata.uns:
    - `guide_library`: Guide library metadata

    Examples
    --------
    >>> import scanpy as sc
    >>> import perturbio as pt
    >>> adata = sc.read_h5ad("cropseq_data.h5ad")
    >>> pt.guides.extract(adata, guide_file="guides.csv")
    >>> print(adata.obs['perturbation'].value_counts())
    """
    adata = adata.copy() if copy else adata

    # Load guide library
    if guide_library is None:
        if guide_file is None:
            raise ValueError(
                "Must provide either guide_file or guide_library\n"
                "Example:\n"
                "  pt.guides.extract(adata, guide_file='guides.csv')"
            )
        guide_library = load_guide_library(guide_file)

    # Store guide library in adata
    adata.uns['guide_library'] = guide_library

    # Extract guides from cells
    print("\n[2/5] Extracting CRISPR guide barcodes...")
    print(f"  ✓ Loaded guide library: {len(guide_library)} guides")

    # Assign guides to cells
    guide_assignments = assign_guides_to_cells(
        adata,
        guide_library,
        min_umis=min_umis,
        multiplet_threshold=multiplet_threshold,
    )

    # Add results to adata.obs
    adata.obs['guide_identity'] = guide_assignments['guide_identity']
    adata.obs[key_added] = guide_assignments['perturbation']
    adata.obs['guide_umi_count'] = guide_assignments['guide_umi_count']
    adata.obs['guide_confidence'] = guide_assignments['confidence']

    # Print summary
    _print_guide_summary(adata, key_added)

    return adata if copy else None


def assign_guides_to_cells(
    adata: ad.AnnData,
    guide_library: pd.DataFrame,
    min_umis: int = 3,
    multiplet_threshold: float = 0.2,
) -> pd.DataFrame:
    """
    Assign guides to cells using winner-takes-all approach.

    Algorithm:
    1. Find guide RNAs in gene expression matrix
    2. For each cell, get UMI counts for all guides
    3. Assign guide with highest count (if above min_umis)
    4. Flag multiplets (cells with multiple high-count guides)

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    guide_library : DataFrame
        Guide library with columns: guide_id, target_gene
    min_umis : int, default 3
        Minimum UMIs required for assignment
    multiplet_threshold : float, default 0.2
        Fraction of max count for multiplet detection

    Returns
    -------
    DataFrame
        Guide assignments with columns:
        - guide_identity: Assigned guide ID
        - perturbation: Target gene
        - guide_umi_count: UMI count
        - confidence: Assignment confidence level
    """
    # Find which guides are present in the data
    guide_ids = guide_library['guide_id'].values
    guides_in_data = [g for g in guide_ids if g in adata.var_names]

    if len(guides_in_data) == 0:
        # No guides found - provide helpful error
        raise ValueError(
            "No guides detected in dataset\n\n"
            f"Your guide library contains:\n"
            f"  • {', '.join(guide_ids[:5])}{'...' if len(guide_ids) > 5 else ''}\n\n"
            f"Your dataset contains genes like:\n"
            f"  • {', '.join(adata.var_names[:5].tolist())}...\n\n"
            "Possible issues:\n"
            "  ✗ Guide names don't match gene IDs in dataset\n"
            "  ✗ Guides may be in a separate assay\n"
            "  ✗ Wrong guide library file\n\n"
            "Solutions:\n"
            "  • Ensure guide IDs in CSV match gene names in h5ad\n"
            "  • Check if guides are in adata.layers or adata.obsm"
        )

    print(f"  ✓ Found {len(guides_in_data)}/{len(guide_ids)} guides in dataset")

    # Get guide expression matrix (cells × guides)
    guide_indices = [list(adata.var_names).index(g) for g in guides_in_data]
    guide_expr = adata.X[:, guide_indices]

    # Convert to dense if sparse
    if hasattr(guide_expr, 'toarray'):
        guide_expr = guide_expr.toarray()

    # Create mapping from guide_id to target_gene
    guide_to_target = dict(zip(guide_library['guide_id'], guide_library['target_gene']))

    # Initialize results
    n_cells = adata.n_obs
    results = pd.DataFrame({
        'guide_identity': ['unassigned'] * n_cells,
        'perturbation': ['unassigned'] * n_cells,
        'guide_umi_count': [0.0] * n_cells,
        'confidence': ['unassigned'] * n_cells,
    }, index=adata.obs_names)

    # Assign guides cell by cell
    print("  ↻ Assigning guides to cells...")
    for i in tqdm(range(n_cells), desc="Assigning", ncols=70):
        cell_guides = guide_expr[i, :]

        # Get guides with non-zero counts
        nonzero_idx = np.where(cell_guides > 0)[0]

        if len(nonzero_idx) == 0:
            # No guides detected
            continue

        # Get counts and guide names
        counts = cell_guides[nonzero_idx]
        guide_names = [guides_in_data[idx] for idx in nonzero_idx]

        # Find guide with max count
        max_idx = np.argmax(counts)
        max_count = counts[max_idx]
        assigned_guide = guide_names[max_idx]

        # Check for multiplets
        if len(counts) > 1:
            # Find guides with counts close to max
            high_counts = counts >= (max_count * (1 - multiplet_threshold))
            if np.sum(high_counts) > 1:
                # Multiplet detected
                results.iloc[i] = {
                    'guide_identity': 'multiplet',
                    'perturbation': 'multiplet',
                    'guide_umi_count': max_count,
                    'confidence': 'multiplet',
                }
                continue

        # Assign guide
        target_gene = guide_to_target.get(assigned_guide, 'unknown')

        # Determine confidence
        if max_count < min_umis:
            confidence = 'low'
        else:
            confidence = 'high'

        results.iloc[i] = {
            'guide_identity': assigned_guide,
            'perturbation': target_gene,
            'guide_umi_count': max_count,
            'confidence': confidence,
        }

    return results


def detect_multiplets(
    adata: ad.AnnData,
    guide_key: str = 'guide_identity',
    threshold: float = 0.2,
) -> pd.Series:
    """
    Detect cells with multiple guides (multiplets).

    Parameters
    ----------
    adata : AnnData
        Annotated data with guide assignments
    guide_key : str, default 'guide_identity'
        Key in adata.obs containing guide assignments
    threshold : float, default 0.2
        Threshold for multiplet detection

    Returns
    -------
    Series
        Boolean series indicating multiplet cells
    """
    if guide_key not in adata.obs.columns:
        raise ValueError(
            f"Guide key '{guide_key}' not found in adata.obs\n"
            "Run pt.guides.extract() first"
        )

    is_multiplet = adata.obs[guide_key] == 'multiplet'
    return is_multiplet


def _print_guide_summary(adata: ad.AnnData, key: str = 'perturbation') -> None:
    """Print summary of guide assignments."""
    perturbations = adata.obs[key]

    # Count categories
    n_assigned = (perturbations != 'unassigned') & (perturbations != 'multiplet')
    n_assigned = n_assigned.sum()

    n_control = (perturbations.str.contains('control|non-targeting|NTC', case=False, na=False)).sum()
    n_unassigned = (perturbations == 'unassigned').sum()
    n_multiplet = (perturbations == 'multiplet').sum()

    # Count low confidence
    if 'guide_confidence' in adata.obs.columns:
        n_low_conf = (adata.obs['guide_confidence'] == 'low').sum()
    else:
        n_low_conf = 0

    # Print table
    print("\n  Guide Assignment Summary:")
    print("  ┌──────────────────────────┬────────┐")
    print("  │ Category                 │ Count  │")
    print("  ├──────────────────────────┼────────┤")
    print(f"  │ Cells with guides        │ {n_assigned:>6} │")
    if n_control > 0:
        print(f"  │ Non-targeting controls   │ {n_control:>6} │")
    print(f"  │ Unassigned cells         │ {n_unassigned:>6} │")
    if n_multiplet > 0:
        print(f"  │ Multiplets               │ {n_multiplet:>6} │")
    if n_low_conf > 0:
        print(f"  │ Low confidence           │ {n_low_conf:>6} │")
    print("  └──────────────────────────┴────────┘")

    # Warning for low assignment rate
    pct_assigned = 100 * n_assigned / adata.n_obs
    if pct_assigned < 50:
        print(f"\n  ⚠ Warning: Only {pct_assigned:.1f}% of cells have guide assignments")
        print("  → Consider adjusting --min-umis threshold")
        print("  → Check that guide names match gene IDs in dataset")
