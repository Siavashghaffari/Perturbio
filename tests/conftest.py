"""
Pytest fixtures for testing Perturbio.
"""

import pytest
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix


@pytest.fixture
def mock_adata():
    """Create a mock AnnData object for testing."""
    np.random.seed(42)

    # Create fake expression data
    n_cells = 100
    n_genes = 50

    # Sparse matrix with counts
    X = csr_matrix(np.random.poisson(2, size=(n_cells, n_genes)).astype(np.float32))

    # Create gene names
    var_names = [f"Gene_{i}" for i in range(n_genes)]

    # Add some guide RNA genes
    guide_names = [
        "BRCA1_guide1",
        "BRCA1_guide2",
        "MYC_guide1",
        "TP53_guide1",
        "non-targeting_1",
    ]
    var_names[:len(guide_names)] = guide_names

    # Create AnnData
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=[f"Cell_{i}" for i in range(n_cells)]),
        var=pd.DataFrame(index=var_names),
    )

    # Add guide expression (make some cells express each guide)
    for i, guide in enumerate(guide_names):
        guide_idx = list(adata.var_names).index(guide)
        # Assign guide to specific cells
        start = i * 15
        end = start + 15
        if end <= n_cells:
            adata.X[start:end, guide_idx] = np.random.poisson(10, size=end-start)

    return adata


@pytest.fixture
def mock_guide_library():
    """Create a mock guide library DataFrame."""
    guides = pd.DataFrame({
        'guide_id': [
            'BRCA1_guide1',
            'BRCA1_guide2',
            'MYC_guide1',
            'TP53_guide1',
            'non-targeting_1',
        ],
        'target_gene': [
            'BRCA1',
            'BRCA1',
            'MYC',
            'TP53',
            'control',
        ],
        'guide_sequence': [
            'GCACTCAGGAAACAGCTATG',
            'CTGAAGACTGCTCAGTGTAG',
            'GTACTTGGTGAGGCCAGCGC',
            'TACAGCGTGGTGGTGCCTAT',
            'GTAGCGAACGTGTCCGGCGT',
        ],
    })
    return guides


@pytest.fixture
def mock_de_results():
    """Create mock differential expression results."""
    n_genes = 50
    perturbations = ['BRCA1', 'MYC', 'TP53']

    results = []
    for pert in perturbations:
        for i in range(n_genes):
            results.append({
                'perturbation': pert,
                'gene': f'Gene_{i}',
                'log_fc': np.random.normal(-1, 1),
                'pval': np.random.uniform(0, 1),
                'pval_adj': np.random.uniform(0, 1),
                'rank': i + 1,
            })

    return pd.DataFrame(results)
