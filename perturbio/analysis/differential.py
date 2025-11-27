"""
Differential expression analysis for Crop-Seq perturbations.
"""

from typing import Optional, Union, List
import warnings

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from tqdm import tqdm


def differential_expression(
    adata: ad.AnnData,
    groupby: str = 'perturbation',
    control: Union[str, List[str]] = 'non-targeting',
    method: str = 'wilcoxon',
    min_cells: int = 10,
    fdr_threshold: float = 0.05,
    key_added: str = 'perturbio_de',
    copy: bool = False,
) -> Optional[ad.AnnData]:
    """
    Perform differential expression analysis for perturbations vs control.

    This function compares gene expression in perturbed cells versus control cells
    using Wilcoxon rank-sum test (or other methods). It tests each perturbation
    separately against the control group.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with perturbation annotations
    groupby : str, default 'perturbation'
        Column in adata.obs containing perturbation labels
    control : str or list, default 'non-targeting'
        Label(s) for control cells. Can be a string or list of strings.
        Also matches 'control', 'NTC' case-insensitively.
    method : str, default 'wilcoxon'
        Statistical test to use. Options: 'wilcoxon', 't-test', 'logreg'
    min_cells : int, default 10
        Minimum number of cells required for a perturbation to be tested
    fdr_threshold : float, default 0.05
        FDR threshold for significance
    key_added : str, default 'perturbio_de'
        Key for storing results in adata.uns
    copy : bool, default False
        Whether to return a copy of adata or modify in place

    Returns
    -------
    AnnData or None
        If copy=True, returns modified copy of adata.
        If copy=False, modifies adata in place and returns None.

    Notes
    -----
    Results are stored in adata.uns[key_added] as a DataFrame with columns:
    - perturbation: Perturbation label
    - gene: Gene name
    - log_fc: Log fold-change vs control
    - pval: P-value
    - pval_adj: Adjusted p-value (FDR)
    - rank: Rank by significance

    Examples
    --------
    >>> import scanpy as sc
    >>> import perturbio as pt
    >>> adata = sc.read_h5ad("cropseq_data.h5ad")
    >>> pt.guides.extract(adata, guide_file="guides.csv")
    >>> pt.tl.differential_expression(adata, groupby='perturbation')
    >>> de_results = adata.uns['perturbio_de']
    """
    adata = adata.copy() if copy else adata

    # Check that groupby column exists
    if groupby not in adata.obs.columns:
        raise ValueError(
            f"Column '{groupby}' not found in adata.obs\n"
            f"Available columns: {list(adata.obs.columns)}\n"
            "Run pt.guides.extract() first to add perturbation annotations"
        )

    print("\n[3/5] Running differential expression...")

    # Identify control cells
    control_mask = _identify_control_cells(adata, groupby, control)
    n_control = control_mask.sum()

    if n_control == 0:
        raise ValueError(
            f"No control cells found with label '{control}'\n"
            f"Available labels: {adata.obs[groupby].unique()}\n"
            "Specify control label with: control='your_control_label'"
        )

    print(f"  ✓ Identified {n_control:,} control cells")

    # Get perturbations to test
    perturbations = adata.obs[groupby].value_counts()

    # Exclude control labels
    if isinstance(control, str):
        control_labels = [control]
    else:
        control_labels = control
    # Also exclude common control patterns
    for pert in perturbations.index:
        if any(c.lower() in str(pert).lower() for c in ['control', 'non-targeting', 'ntc']):
            if pert not in control_labels:
                control_labels.append(pert)

    # Exclude unassigned and multiplets
    exclude_labels = control_labels + ['unassigned', 'multiplet']
    perturbations_to_test = [
        p for p in perturbations.index
        if p not in exclude_labels and perturbations[p] >= min_cells
    ]

    # Warn about skipped perturbations
    skipped = [p for p in perturbations.index
               if p not in exclude_labels and p not in perturbations_to_test]
    if skipped:
        print(f"  ⚠ Warning: Skipping {len(skipped)} perturbations with <{min_cells} cells:")
        for p in skipped[:5]:  # Show first 5
            print(f"    → {p} ({perturbations[p]} cells)")
        if len(skipped) > 5:
            print(f"    ... and {len(skipped) - 5} more")

    n_test = len(perturbations_to_test)
    if n_test == 0:
        raise ValueError(
            f"No perturbations to test (all have <{min_cells} cells)\n"
            f"Consider lowering min_cells parameter"
        )

    print(f"  ↻ Testing {n_test} perturbations vs control...")

    # Run differential expression for all perturbations
    de_results = rank_genes_per_perturbation(
        adata,
        perturbations=perturbations_to_test,
        control_mask=control_mask,
        groupby=groupby,
        method=method,
    )

    # Store results
    adata.uns[key_added] = de_results

    # Print summary
    _print_de_summary(de_results, fdr_threshold)

    return adata if copy else None


def rank_genes_per_perturbation(
    adata: ad.AnnData,
    perturbations: List[str],
    control_mask: pd.Series,
    groupby: str = 'perturbation',
    method: str = 'wilcoxon',
) -> pd.DataFrame:
    """
    Rank genes for each perturbation vs control using scanpy.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    perturbations : list
        List of perturbation labels to test
    control_mask : Series
        Boolean mask for control cells
    groupby : str
        Column name for grouping
    method : str
        Statistical test method

    Returns
    -------
    DataFrame
        Differential expression results for all perturbations
    """
    all_results = []

    # Test each perturbation separately
    for pert in tqdm(perturbations, desc="Testing", ncols=70):
        # Create subset with this perturbation + controls
        pert_mask = adata.obs[groupby] == pert
        subset_mask = pert_mask | control_mask
        adata_subset = adata[subset_mask, :].copy()

        # Create binary label for DE test
        adata_subset.obs['_de_group'] = 'control'
        adata_subset.obs.loc[adata_subset.obs[groupby] == pert, '_de_group'] = 'perturbation'

        # Run DE test using scanpy
        try:
            sc.tl.rank_genes_groups(
                adata_subset,
                groupby='_de_group',
                groups=['perturbation'],
                reference='control',
                method=method,
                use_raw=False,
                key_added='_de_temp',
            )

            # Extract results
            result = sc.get.rank_genes_groups_df(
                adata_subset,
                group='perturbation',
                key='_de_temp',
            )

            # Add perturbation label
            result.insert(0, 'perturbation', pert)

            # Rename columns for consistency
            result = result.rename(columns={
                'names': 'gene',
                'logfoldchanges': 'log_fc',
                'pvals': 'pval',
                'pvals_adj': 'pval_adj',
            })

            # Add rank
            result['rank'] = range(1, len(result) + 1)

            all_results.append(result)

        except Exception as e:
            warnings.warn(f"Failed to test perturbation '{pert}': {str(e)}", UserWarning)
            continue

    # Combine all results
    if len(all_results) == 0:
        raise ValueError("No perturbations could be tested successfully")

    combined = pd.concat(all_results, ignore_index=True)

    # Select and order columns
    columns = ['perturbation', 'gene', 'log_fc', 'pval', 'pval_adj', 'rank']
    combined = combined[columns]

    return combined


def _identify_control_cells(
    adata: ad.AnnData,
    groupby: str,
    control: Union[str, List[str]],
) -> pd.Series:
    """Identify control cells based on label(s)."""
    if isinstance(control, str):
        control_labels = [control]
    else:
        control_labels = control

    # Create mask for control cells
    control_mask = pd.Series(False, index=adata.obs_names)

    for label in control_labels:
        # Exact match
        control_mask |= (adata.obs[groupby] == label)

        # Case-insensitive partial match for common patterns
        for pattern in [label, 'control', 'non-targeting', 'ntc']:
            pattern_mask = adata.obs[groupby].str.contains(
                pattern, case=False, na=False, regex=False
            )
            control_mask |= pattern_mask

    return control_mask


def _print_de_summary(de_results: pd.DataFrame, fdr_threshold: float = 0.05) -> None:
    """Print summary of differential expression results."""
    n_sig = (de_results['pval_adj'] < fdr_threshold).sum()
    print(f"\n  ✓ Found {n_sig:,} significantly perturbed genes (FDR < {fdr_threshold})")

    # Get top hits for each perturbation
    top_hits = []
    for pert in de_results['perturbation'].unique():
        pert_data = de_results[de_results['perturbation'] == pert]
        sig_genes = (pert_data['pval_adj'] < fdr_threshold).sum()

        if len(pert_data) > 0:
            top_gene = pert_data.iloc[0]
            top_hits.append({
                'perturbation': pert,
                'cells': 'N/A',  # Will be filled in by caller if available
                'sig_genes': sig_genes,
                'top_gene': top_gene['gene'],
                'log_fc': top_gene['log_fc'],
            })

    if len(top_hits) > 0:
        # Show top 5 perturbations by number of significant genes
        top_hits_df = pd.DataFrame(top_hits).sort_values('sig_genes', ascending=False)

        print("\n  Perturbation Results:")
        print("  ┌─────────────────┬───────┬──────────┬─────────────┐")
        print("  │ Perturbation    │ Cells │ Sig Gene │ Top Gene    │")
        print("  ├─────────────────┼───────┼──────────┼─────────────┤")

        for _, row in top_hits_df.head(10).iterrows():
            pert_short = row['perturbation'][:15]
            gene_with_fc = f"{row['top_gene']} ↓{abs(row['log_fc']):.1f}"
            print(f"  │ {pert_short:<15} │  N/A  │ {row['sig_genes']:>6}   │ {gene_with_fc:<11} │")

        if len(top_hits_df) > 10:
            print("  │ ...             │  ...  │   ...    │ ...         │")

        print("  └─────────────────┴───────┴──────────┴─────────────┘")
