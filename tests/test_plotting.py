"""
Tests for plotting functionality.
"""

import pytest
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from perturbio.plotting.core import volcano, umap, perturbation_counts


def test_volcano_with_adata(mock_adata, mock_guide_library, tmp_path):
    """Test volcano plot with AnnData object."""
    from perturbio.guides.extraction import extract
    from perturbio.analysis.differential import differential_expression

    # Prepare data
    guide_file = tmp_path / "guides.csv"
    mock_guide_library.to_csv(guide_file, index=False)
    extract(mock_adata, guide_file=guide_file, min_umis=1)
    differential_expression(mock_adata, control='control', min_cells=5)

    # Create volcano plot
    fig = volcano(adata=mock_adata, perturbation='BRCA1')

    assert fig is not None
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_volcano_with_de_results(mock_de_results):
    """Test volcano plot with DE results DataFrame."""
    fig = volcano(de_results=mock_de_results, perturbation='BRCA1')

    assert fig is not None
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_volcano_custom_params(mock_de_results):
    """Test volcano plot with custom parameters."""
    fig = volcano(
        de_results=mock_de_results,
        perturbation='BRCA1',
        fdr_threshold=0.01,
        log_fc_threshold=1.0,
        label_top=5,
        figsize=(10, 8)
    )

    assert fig is not None
    plt.close(fig)


def test_volcano_save_file(mock_de_results, tmp_path):
    """Test volcano plot saving to file."""
    output_file = tmp_path / "volcano.png"

    fig = volcano(
        de_results=mock_de_results,
        perturbation='BRCA1',
        save=str(output_file)
    )

    assert output_file.exists()
    plt.close(fig)


@pytest.mark.skip(reason="UMAP causes TensorFlow import crash in test environment")
def test_umap_basic(mock_adata, mock_guide_library, tmp_path):
    """Test basic UMAP plot."""
    from perturbio.guides.extraction import extract
    import scanpy as sc

    # Prepare data
    guide_file = tmp_path / "guides.csv"
    mock_guide_library.to_csv(guide_file, index=False)
    extract(mock_adata, guide_file=guide_file, min_umis=1)

    # Add UMAP coordinates
    sc.pp.pca(mock_adata, n_comps=5)
    sc.pp.neighbors(mock_adata)
    sc.tl.umap(mock_adata)

    # Create UMAP plot
    fig = umap(mock_adata, color='perturbation')

    assert fig is not None
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


@pytest.mark.skip(reason="UMAP causes TensorFlow import crash in test environment")
def test_umap_save_file(mock_adata, mock_guide_library, tmp_path):
    """Test UMAP plot saving."""
    from perturbio.guides.extraction import extract
    import scanpy as sc

    guide_file = tmp_path / "guides.csv"
    mock_guide_library.to_csv(guide_file, index=False)
    extract(mock_adata, guide_file=guide_file, min_umis=1)

    sc.pp.pca(mock_adata, n_comps=5)
    sc.pp.neighbors(mock_adata)
    sc.tl.umap(mock_adata)

    output_file = tmp_path / "umap.png"
    fig = umap(mock_adata, color='perturbation', save=str(output_file))

    assert output_file.exists()
    plt.close(fig)


def test_perturbation_counts_basic(mock_adata, mock_guide_library, tmp_path):
    """Test basic perturbation counts plot."""
    from perturbio.guides.extraction import extract

    guide_file = tmp_path / "guides.csv"
    mock_guide_library.to_csv(guide_file, index=False)
    extract(mock_adata, guide_file=guide_file, min_umis=1)

    fig = perturbation_counts(mock_adata, groupby='perturbation')

    assert fig is not None
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_perturbation_counts_min_cells(mock_adata, mock_guide_library, tmp_path):
    """Test perturbation counts with min_cells filter."""
    from perturbio.guides.extraction import extract

    guide_file = tmp_path / "guides.csv"
    mock_guide_library.to_csv(guide_file, index=False)
    extract(mock_adata, guide_file=guide_file, min_umis=1)

    fig = perturbation_counts(mock_adata, groupby='perturbation', min_cells=10)

    assert fig is not None
    plt.close(fig)


def test_perturbation_counts_save(mock_adata, mock_guide_library, tmp_path):
    """Test perturbation counts saving."""
    from perturbio.guides.extraction import extract

    guide_file = tmp_path / "guides.csv"
    mock_guide_library.to_csv(guide_file, index=False)
    extract(mock_adata, guide_file=guide_file, min_umis=1)

    output_file = tmp_path / "counts.png"
    fig = perturbation_counts(mock_adata, save=str(output_file))

    assert output_file.exists()
    plt.close(fig)


def test_volcano_no_data():
    """Test volcano plot with no data raises error."""
    with pytest.raises(ValueError, match="Must provide either"):
        volcano(perturbation='BRCA1')


def test_volcano_no_perturbation(mock_de_results):
    """Test volcano plot without perturbation specified."""
    with pytest.raises((ValueError, KeyError)):
        volcano(de_results=mock_de_results, perturbation='NONEXISTENT')
