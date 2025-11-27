"""
Tests for differential expression analysis.
"""

import pytest
import pandas as pd
from perturbio.analysis.differential import differential_expression


def test_differential_expression(mock_adata, mock_guide_library, tmp_path):
    """Test differential expression analysis."""
    # First extract guides
    from perturbio.guides.extraction import extract

    guide_file = tmp_path / "guides.csv"
    mock_guide_library.to_csv(guide_file, index=False)
    extract(mock_adata, guide_file=guide_file, min_umis=1)

    # Run differential expression
    differential_expression(
        mock_adata,
        groupby='perturbation',
        control='control',
        min_cells=5,
    )

    # Check that results were added
    assert 'perturbio_de' in mock_adata.uns

    # Check results format
    de_results = mock_adata.uns['perturbio_de']
    assert isinstance(de_results, pd.DataFrame)

    # Check required columns
    required_cols = ['perturbation', 'gene', 'log_fc', 'pval', 'pval_adj', 'rank']
    for col in required_cols:
        assert col in de_results.columns


def test_differential_expression_no_guides(mock_adata):
    """Test that DE raises error without guide extraction."""
    with pytest.raises(ValueError, match="not found in adata.obs"):
        differential_expression(mock_adata, groupby='perturbation')


def test_differential_expression_no_control(mock_adata, mock_guide_library, tmp_path):
    """Test error when no control cells found."""
    from perturbio.guides.extraction import extract

    guide_file = tmp_path / "guides.csv"
    mock_guide_library.to_csv(guide_file, index=False)
    extract(mock_adata, guide_file=guide_file)

    with pytest.raises(ValueError, match="No control cells found"):
        differential_expression(
            mock_adata,
            groupby='perturbation',
            control='NONEXISTENT',
        )
