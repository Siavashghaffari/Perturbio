"""
Tests for guide extraction module.
"""

import pytest
import pandas as pd
from perturbio.guides.extraction import extract, assign_guides_to_cells


def test_extract_guides(mock_adata, mock_guide_library, tmp_path):
    """Test guide extraction from AnnData."""
    # Save guide library to temp file
    guide_file = tmp_path / "guides.csv"
    mock_guide_library.to_csv(guide_file, index=False)

    # Extract guides
    extract(mock_adata, guide_file=guide_file, min_umis=3)

    # Check that annotations were added
    assert 'guide_identity' in mock_adata.obs.columns
    assert 'perturbation' in mock_adata.obs.columns
    assert 'guide_umi_count' in mock_adata.obs.columns
    assert 'guide_confidence' in mock_adata.obs.columns

    # Check that some cells were assigned
    assigned = mock_adata.obs['perturbation'] != 'unassigned'
    assert assigned.sum() > 0

    # Check that guide library was stored
    assert 'guide_library' in mock_adata.uns


def test_assign_guides_to_cells(mock_adata, mock_guide_library):
    """Test guide assignment algorithm."""
    assignments = assign_guides_to_cells(
        mock_adata,
        mock_guide_library,
        min_umis=3,
    )

    # Check that we got a DataFrame back
    assert isinstance(assignments, pd.DataFrame)
    assert len(assignments) == mock_adata.n_obs

    # Check required columns
    required_cols = ['guide_identity', 'perturbation', 'guide_umi_count', 'confidence']
    for col in required_cols:
        assert col in assignments.columns


def test_extract_guides_no_library(mock_adata):
    """Test that extract raises error without guide library."""
    with pytest.raises(ValueError, match="Must provide either guide_file or guide_library"):
        extract(mock_adata)


def test_extract_guides_no_matches(mock_adata):
    """Test error when no guides match the data."""
    # Create guide library with non-existent guides
    fake_guides = pd.DataFrame({
        'guide_id': ['FAKE_guide1', 'FAKE_guide2'],
        'target_gene': ['FAKE1', 'FAKE2'],
    })

    with pytest.raises(ValueError, match="No guides detected in dataset"):
        assign_guides_to_cells(mock_adata, fake_guides)
