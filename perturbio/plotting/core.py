"""
Core plotting functions for Crop-Seq analysis.
"""

from typing import Optional, Union
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import anndata as ad


def volcano(
    adata: Optional[ad.AnnData] = None,
    de_results: Optional[pd.DataFrame] = None,
    perturbation: str = None,
    fdr_threshold: float = 0.05,
    log_fc_threshold: float = 0.5,
    label_top: int = 10,
    figsize: tuple = (8, 6),
    save: Optional[Union[str, Path]] = None,
) -> plt.Figure:
    """
    Create a volcano plot for differential expression results.

    Parameters
    ----------
    adata : AnnData, optional
        Annotated data with DE results in adata.uns['perturbio_de']
    de_results : DataFrame, optional
        Differential expression results DataFrame
    perturbation : str
        Perturbation to plot
    fdr_threshold : float, default 0.05
        FDR threshold for significance
    log_fc_threshold : float, default 0.5
        Log fold-change threshold for highlighting
    label_top : int, default 10
        Number of top genes to label
    figsize : tuple, default (8, 6)
        Figure size
    save : str or Path, optional
        Path to save figure

    Returns
    -------
    Figure
        Matplotlib figure object

    Examples
    --------
    >>> import perturbio as pt
    >>> fig = pt.pl.volcano(adata, perturbation='BRCA1_guide1')
    >>> plt.show()
    """
    # Get DE results
    if de_results is None:
        if adata is None:
            raise ValueError("Must provide either adata or de_results")
        if 'perturbio_de' not in adata.uns:
            raise ValueError(
                "No DE results found in adata.uns['perturbio_de']\n"
                "Run pt.tl.differential_expression() first"
            )
        de_results = adata.uns['perturbio_de']

    # Filter for perturbation
    if perturbation is None:
        raise ValueError("Must specify perturbation to plot")

    data = de_results[de_results['perturbation'] == perturbation].copy()

    if len(data) == 0:
        raise ValueError(f"No results found for perturbation '{perturbation}'")

    # Calculate -log10(pval_adj)
    data['-log10_pval'] = -np.log10(data['pval_adj'].clip(lower=1e-300))

    # Determine significance
    data['significant'] = (
        (data['pval_adj'] < fdr_threshold) &
        (np.abs(data['log_fc']) > log_fc_threshold)
    )

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Plot non-significant points
    nonsig = data[~data['significant']]
    ax.scatter(
        nonsig['log_fc'],
        nonsig['-log10_pval'],
        c='lightgray',
        s=10,
        alpha=0.5,
        label='Not significant',
    )

    # Plot significant points
    sig = data[data['significant']]
    if len(sig) > 0:
        ax.scatter(
            sig['log_fc'],
            sig['-log10_pval'],
            c='red',
            s=20,
            alpha=0.7,
            label=f'Significant (FDR < {fdr_threshold})',
        )

        # Label top genes
        top_genes = sig.nlargest(label_top, '-log10_pval')
        for _, gene in top_genes.iterrows():
            ax.annotate(
                gene['gene'],
                xy=(gene['log_fc'], gene['-log10_pval']),
                xytext=(5, 5),
                textcoords='offset points',
                fontsize=8,
                alpha=0.8,
            )

    # Add threshold lines
    ax.axhline(
        y=-np.log10(fdr_threshold),
        color='gray',
        linestyle='--',
        linewidth=1,
        alpha=0.5,
        label=f'FDR = {fdr_threshold}',
    )
    ax.axvline(x=0, color='gray', linestyle='-', linewidth=0.5, alpha=0.3)

    # Labels and title
    ax.set_xlabel('Log2 Fold Change', fontsize=12)
    ax.set_ylabel('-Log10 Adjusted P-value', fontsize=12)
    ax.set_title(f'Volcano Plot: {perturbation}', fontsize=14, fontweight='bold')
    ax.legend(loc='upper right', frameon=True, fontsize=9)

    # Style
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(True, alpha=0.2)

    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=300, bbox_inches='tight')
        print(f"✓ Saved volcano plot: {save}")

    return fig


def umap(
    adata: ad.AnnData,
    color: str = 'perturbation',
    legend_loc: str = 'right margin',
    figsize: tuple = (8, 6),
    save: Optional[Union[str, Path]] = None,
) -> plt.Figure:
    """
    Create UMAP plot colored by perturbations.

    Parameters
    ----------
    adata : AnnData
        Annotated data with UMAP coordinates
    color : str, default 'perturbation'
        Column in adata.obs to color by
    legend_loc : str, default 'right margin'
        Legend location
    figsize : tuple, default (8, 6)
        Figure size
    save : str or Path, optional
        Path to save figure

    Returns
    -------
    Figure
        Matplotlib figure object

    Examples
    --------
    >>> import scanpy as sc
    >>> import perturbio as pt
    >>> sc.tl.umap(adata)
    >>> fig = pt.pl.umap(adata, color='perturbation')
    """
    # Check for UMAP coordinates
    if 'X_umap' not in adata.obsm:
        raise ValueError(
            "No UMAP coordinates found in adata.obsm['X_umap']\n"
            "Run sc.tl.umap(adata) first"
        )

    # Check color column
    if color not in adata.obs.columns:
        raise ValueError(
            f"Column '{color}' not found in adata.obs\n"
            f"Available: {list(adata.obs.columns)}"
        )

    # Get UMAP coordinates
    umap_coords = adata.obsm['X_umap']
    colors_array = adata.obs[color]

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Get unique categories
    categories = colors_array.unique()
    n_cats = len(categories)

    # Use a good color palette
    if n_cats <= 10:
        palette = sns.color_palette("tab10", n_cats)
    elif n_cats <= 20:
        palette = sns.color_palette("tab20", n_cats)
    else:
        palette = sns.color_palette("husl", n_cats)

    # Plot each category
    for i, cat in enumerate(categories):
        mask = colors_array == cat
        ax.scatter(
            umap_coords[mask, 0],
            umap_coords[mask, 1],
            c=[palette[i]],
            label=cat,
            s=10,
            alpha=0.7,
        )

    ax.set_xlabel('UMAP 1', fontsize=12)
    ax.set_ylabel('UMAP 2', fontsize=12)
    ax.set_title(f'UMAP colored by {color}', fontsize=14, fontweight='bold')

    # Legend
    if legend_loc == 'right margin':
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=True, fontsize=8)
    else:
        ax.legend(loc=legend_loc, frameon=True, fontsize=8)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=300, bbox_inches='tight')
        print(f"✓ Saved UMAP plot: {save}")

    return fig


def perturbation_counts(
    adata: ad.AnnData,
    groupby: str = 'perturbation',
    min_cells: int = 0,
    figsize: tuple = (10, 6),
    save: Optional[Union[str, Path]] = None,
) -> plt.Figure:
    """
    Create bar chart of cell counts per perturbation.

    Parameters
    ----------
    adata : AnnData
        Annotated data with perturbation labels
    groupby : str, default 'perturbation'
        Column in adata.obs to count
    min_cells : int, default 0
        Only show perturbations with at least this many cells
    figsize : tuple, default (10, 6)
        Figure size
    save : str or Path, optional
        Path to save figure

    Returns
    -------
    Figure
        Matplotlib figure object

    Examples
    --------
    >>> import perturbio as pt
    >>> fig = pt.pl.perturbation_counts(adata)
    """
    if groupby not in adata.obs.columns:
        raise ValueError(
            f"Column '{groupby}' not found in adata.obs\n"
            f"Available: {list(adata.obs.columns)}"
        )

    # Count cells per perturbation
    counts = adata.obs[groupby].value_counts()

    # Filter by min_cells
    if min_cells > 0:
        counts = counts[counts >= min_cells]

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Create bar plot
    x_pos = np.arange(len(counts))
    bars = ax.bar(x_pos, counts.values, color='steelblue', alpha=0.7)

    # Highlight control bars
    for i, (label, count) in enumerate(counts.items()):
        if any(c in str(label).lower() for c in ['control', 'non-targeting', 'ntc']):
            bars[i].set_color('lightgray')

    # Labels
    ax.set_xlabel('Perturbation', fontsize=12)
    ax.set_ylabel('Number of Cells', fontsize=12)
    ax.set_title('Cell Counts per Perturbation', fontsize=14, fontweight='bold')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(counts.index, rotation=45, ha='right', fontsize=9)

    # Add count labels on bars
    for i, (label, count) in enumerate(counts.items()):
        ax.text(
            i, count + max(counts) * 0.01,
            str(count),
            ha='center',
            va='bottom',
            fontsize=8,
        )

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', alpha=0.2)

    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=300, bbox_inches='tight')
        print(f"✓ Saved perturbation counts plot: {save}")

    return fig
