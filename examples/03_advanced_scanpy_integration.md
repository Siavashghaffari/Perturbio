# Advanced Scanpy Integration

**Time:** 15-20 minutes
**Level:** Advanced

## Overview

This tutorial shows how to integrate Perturbio seamlessly into your existing scanpy workflows. You'll learn how to use Perturbio's low-level API for maximum flexibility.

## Learning Objectives

- Use Perturbio functions directly with scanpy
- Create custom analysis workflows
- Combine Perturbio with other single-cell tools
- Optimize performance for large datasets

---

## Part 1: Direct Function Calls

Perturbio provides low-level functions that work directly with AnnData objects, following scanpy conventions.

```python
import scanpy as sc
import perturbio as pt
import numpy as np
import pandas as pd

# Load your data
adata = sc.read_h5ad("cropseq_data.h5ad")
```

### Standard Scanpy Preprocessing

```python
# Quality control
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Calculate QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

# Filter low-quality cells
adata = adata[adata.obs['pct_counts_mt'] < 20, :]

# Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

print(f"After QC: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
```

### Extract Guides (In-Place)

```python
# Extract guides - modifies adata in place (scanpy-style)
pt.guides.extract(
    adata,
    guide_file="guides.csv",
    min_umis=3,
    key_added="perturbation"  # Column name in adata.obs
)

# Check what was added
print(adata.obs.columns)
# [..., 'perturbation', 'guide_identity', 'guide_umi_count']
```

### Dimensionality Reduction with Perturbations

```python
# Find variable genes
sc.pp.highly_variable_genes(adata, n_top_genes=2000)

# PCA
sc.tl.pca(adata, n_comps=50)

# Batch correction (if needed)
# import scanpy.external as sce
# sce.pp.harmony_integrate(adata, 'batch')

# Neighbors and UMAP
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
sc.tl.umap(adata, min_dist=0.3)

# Clustering
sc.tl.leiden(adata, resolution=0.5)

# Plot UMAP colored by perturbation
sc.pl.umap(adata, color=['perturbation', 'leiden'], ncols=2)
```

---

## Part 2: Differential Expression Analysis

### Using Perturbio's DE Function

```python
# Run DE analysis (in-place)
pt.tl.differential_expression(
    adata,
    groupby='perturbation',
    control='non-targeting',
    method='wilcoxon',
    min_cells=20,
    fdr_threshold=0.05,
    key_added='perturbio_de'  # Store results here
)

# Access results
de_results = adata.uns['perturbio_de']
print(de_results.head())
```

### Custom DE Analysis Per Cluster

```python
# Analyze perturbation effects within each cluster
for cluster in adata.obs['leiden'].unique():
    print(f"\\nCluster {cluster}:")

    # Subset to this cluster
    adata_cluster = adata[adata.obs['leiden'] == cluster, :]

    # Run DE for perturbations in this cluster
    pt.tl.differential_expression(
        adata_cluster,
        groupby='perturbation',
        control='non-targeting',
        key_added=f'de_cluster_{cluster}'
    )

    # Get results
    cluster_results = adata_cluster.uns[f'de_cluster_{cluster}']
    print(f"  Tested {len(cluster_results['perturbation'].unique())} perturbations")
```

---

## Part 3: Advanced Visualizations

### Custom Volcano Plot

```python
import matplotlib.pyplot as plt

# Create custom volcano plot
def custom_volcano(de_results, perturbation, fdr=0.05, fc_thresh=0.5):
    # Filter for perturbation
    data = de_results[de_results['perturbation'] == perturbation].copy()

    # Calculate -log10(p)
    data['-log10p'] = -np.log10(data['pval_adj'].clip(lower=1e-300))

    # Significance categories
    data['category'] = 'Not significant'
    data.loc[(data['pval_adj'] < fdr) & (data['log_fc'] > fc_thresh), 'category'] = 'Up'
    data.loc[(data['pval_adj'] < fdr) & (data['log_fc'] < -fc_thresh), 'category'] = 'Down'

    # Plot
    fig, ax = plt.subplots(figsize=(10, 8))

    colors = {'Not significant': 'gray', 'Up': 'red', 'Down': 'blue'}
    for cat, color in colors.items():
        subset = data[data['category'] == cat]
        ax.scatter(subset['log_fc'], subset['-log10p'],
                  c=color, label=cat, s=20, alpha=0.6)

    # Add labels for top genes
    top_genes = data.nlargest(10, '-log10p')
    for _, gene in top_genes.iterrows():
        ax.annotate(gene['gene'],
                   xy=(gene['log_fc'], gene['-log10p']),
                   xytext=(5, 5), textcoords='offset points',
                   fontsize=8, alpha=0.8)

    # Styling
    ax.axhline(-np.log10(fdr), ls='--', c='gray', alpha=0.5)
    ax.axvline(fc_thresh, ls='--', c='gray', alpha=0.5)
    ax.axvline(-fc_thresh, ls='--', c='gray', alpha=0.5)

    ax.set_xlabel('Log2 Fold Change', fontsize=12)
    ax.set_ylabel('-Log10 Adjusted P-value', fontsize=12)
    ax.set_title(f'{perturbation} Volcano Plot', fontsize=14, weight='bold')
    ax.legend()
    ax.grid(alpha=0.3)

    return fig

# Use it
fig = custom_volcano(adata.uns['perturbio_de'], 'BRCA1', fdr=0.05)
plt.show()
```

### UMAP with Perturbation Overlays

```python
# Create multi-panel UMAP showing specific perturbations
perturbations_to_show = ['BRCA1', 'MYC', 'TP53', 'non-targeting']

fig, axes = plt.subplots(2, 2, figsize=(14, 12))
axes = axes.flatten()

for i, pert in enumerate(perturbations_to_show):
    # Create binary mask
    mask = adata.obs['perturbation'] == pert

    # Plot background cells
    axes[i].scatter(
        adata.obsm['X_umap'][~mask, 0],
        adata.obsm['X_umap'][~mask, 1],
        c='lightgray', s=5, alpha=0.3
    )

    # Plot cells with this perturbation
    if mask.sum() > 0:
        axes[i].scatter(
            adata.obsm['X_umap'][mask, 0],
            adata.obsm['X_umap'][mask, 1],
            c='red', s=20, alpha=0.8, label=pert
        )

    axes[i].set_title(f'{pert} ({mask.sum()} cells)', fontsize=12, weight='bold')
    axes[i].set_xlabel('UMAP 1')
    axes[i].set_ylabel('UMAP 2')
    axes[i].legend()

plt.tight_layout()
plt.show()
```

---

## Part 4: Performance Optimization

### For Large Datasets (>100K cells)

```python
# Use sparse matrices
adata.X = sc.sparse.csr_matrix(adata.X)

# Process in batches
from tqdm import tqdm

def batch_differential_expression(adata, perturbations, batch_size=10):
    \"\"\"Run DE in batches to save memory.\"\"\"
    all_results = []

    for i in tqdm(range(0, len(perturbations), batch_size)):
        batch = perturbations[i:i+batch_size]

        # Subset to these perturbations + controls
        mask = adata.obs['perturbation'].isin(batch + ['non-targeting'])
        adata_batch = adata[mask, :]

        # Run DE
        pt.tl.differential_expression(
            adata_batch,
            groupby='perturbation',
            control='non-targeting'
        )

        all_results.append(adata_batch.uns['perturbio_de'])

    # Combine results
    return pd.concat(all_results, ignore_index=True)
```

### Parallel Processing

```python
from joblib import Parallel, delayed

def analyze_perturbation(adata, perturbation):
    \"\"\"Analyze single perturbation.\"\"\"
    # Subset
    mask = (adata.obs['perturbation'] == perturbation) | \\
           (adata.obs['perturbation'] == 'non-targeting')
    adata_sub = adata[mask, :].copy()

    # Run DE
    pt.tl.differential_expression(
        adata_sub,
        groupby='perturbation',
        control='non-targeting'
    )

    return adata_sub.uns['perturbio_de']

# Run in parallel
perturbations = adata.obs['perturbation'].unique()
perturbations = [p for p in perturbations if p != 'non-targeting']

results = Parallel(n_jobs=4)(
    delayed(analyze_perturbation)(adata, pert)
    for pert in perturbations
)

# Combine
combined_results = pd.concat(results, ignore_index=True)
```

---

## Part 5: Integration with Other Tools

### Gene Set Enrichment with GSEApy

```python
# After DE analysis, run GSEA
import gseapy as gp

# Get ranked gene list for a perturbation
pert_results = adata.uns['perturbio_de']
pert_results = pert_results[pert_results['perturbation'] == 'BRCA1']

# Create ranked list
ranked_genes = pert_results.set_index('gene')['log_fc'].sort_values(ascending=False)

# Run GSEA
enr = gp.prerank(
    rnk=ranked_genes,
    gene_sets='KEGG_2021_Human',
    outdir='gsea_results',
    verbose=True
)
```

### Trajectory Analysis with PAGA

```python
# Identify trajectory affected by perturbations
sc.tl.paga(adata, groups='leiden')
sc.pl.paga(adata, color='perturbation')

# Compute pseudotime
adata.uns['iroot'] = np.flatnonzero(adata.obs['leiden'] == '0')[0]
sc.tl.dpt(adata)

# Plot pseudotime vs perturbation
sc.pl.scatter(adata, x='dpt_pseudotime', y='perturbation', color='perturbation')
```

---

## Part 6: Export for Other Tools

### Export to Seurat (R)

```python
# Save for Seurat
adata.write('cropseq_with_perturbations.h5ad')

# In R:
# library(Seurat)
# library(SeuratDisk)
# Convert('cropseq_with_perturbations.h5ad', dest='h5seurat')
# seurat_obj <- LoadH5Seurat('cropseq_with_perturbations.h5seurat')
```

### Export Results Tables

```python
# Export DE results to Excel with multiple sheets
with pd.ExcelWriter('perturbio_results.xlsx') as writer:
    # Summary
    summary = adata.uns['perturbio_de'].groupby('perturbation').agg({
        'gene': 'count',
        'pval_adj': lambda x: (x < 0.05).sum()
    })
    summary.columns = ['Total Genes', 'Significant Genes']
    summary.to_excel(writer, sheet_name='Summary')

    # Full results
    adata.uns['perturbio_de'].to_excel(writer, sheet_name='DE_Results', index=False)

    # Top hits per perturbation
    for pert in adata.uns['perturbio_de']['perturbation'].unique()[:10]:
        subset = adata.uns['perturbio_de'][
            adata.uns['perturbio_de']['perturbation'] == pert
        ].head(50)
        subset.to_excel(writer, sheet_name=f'Top_{pert}'[:31], index=False)
```

---

## Summary

You've learned how to:

✅ Integrate Perturbio into scanpy workflows
✅ Use low-level Perturbio functions directly
✅ Create custom visualizations
✅ Optimize for large datasets
✅ Combine with other single-cell tools
✅ Export results for downstream analysis

## Tips and Best Practices

1. **Always normalize your data** before guide extraction
2. **Check guide assignment rates** - aim for >80% assignment
3. **Use enough cells per perturbation** - minimum 20-30 for robust DE
4. **Validate target knockdown** - check if target gene shows reduced expression
5. **Control for batch effects** if analyzing multiple experiments

Happy analyzing! 🧬
