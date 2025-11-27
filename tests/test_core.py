"""
Tests for CropSeqAnalyzer core class.
"""

import pytest
from perturbio.core import CropSeqAnalyzer


def test_analyzer_init_with_adata(mock_adata):
    """Test initializing analyzer with AnnData object."""
    analyzer = CropSeqAnalyzer(mock_adata)

    assert analyzer.adata is mock_adata
    assert analyzer.results is not None
    assert not analyzer._guides_extracted
    assert not analyzer._de_run


def test_analyzer_extract_guides(mock_adata, mock_guide_library, tmp_path):
    """Test guide extraction via analyzer."""
    guide_file = tmp_path / "guides.csv"
    mock_guide_library.to_csv(guide_file, index=False)

    analyzer = CropSeqAnalyzer(mock_adata)
    analyzer.extract_guides(guide_file=guide_file)

    assert analyzer._guides_extracted
    assert 'perturbation' in analyzer.adata.obs.columns
    assert analyzer.results.guide_assignments is not None


def test_analyzer_differential_expression(mock_adata, mock_guide_library, tmp_path):
    """Test differential expression via analyzer."""
    guide_file = tmp_path / "guides.csv"
    mock_guide_library.to_csv(guide_file, index=False)

    analyzer = CropSeqAnalyzer(mock_adata)
    analyzer.extract_guides(guide_file=guide_file, min_umis=1)
    analyzer.differential_expression(control_label='control', min_cells=5)

    assert analyzer._de_run
    assert analyzer.results.differential_expression is not None


def test_analyzer_run(mock_adata, mock_guide_library, tmp_path):
    """Test complete analysis pipeline."""
    guide_file = tmp_path / "guides.csv"
    mock_guide_library.to_csv(guide_file, index=False)

    analyzer = CropSeqAnalyzer(mock_adata)
    results = analyzer.run(guide_file=guide_file, min_cells=5)

    assert results is not None
    assert analyzer._guides_extracted
    assert analyzer._de_run


def test_analyzer_top_hits(mock_adata, mock_guide_library, tmp_path):
    """Test getting top hits."""
    guide_file = tmp_path / "guides.csv"
    mock_guide_library.to_csv(guide_file, index=False)

    analyzer = CropSeqAnalyzer(mock_adata)
    analyzer.run(guide_file=guide_file, min_cells=5)

    # Get perturbations
    perturbations = analyzer.results.differential_expression.perturbations

    if len(perturbations) > 0:
        top_hits = analyzer.top_hits(perturbations[0], n=10)
        assert len(top_hits) <= 10
        assert 'gene' in top_hits.columns


def test_analyzer_export(mock_adata, mock_guide_library, tmp_path):
    """Test exporting results."""
    guide_file = tmp_path / "guides.csv"
    mock_guide_library.to_csv(guide_file, index=False)

    analyzer = CropSeqAnalyzer(mock_adata)
    analyzer.run(guide_file=guide_file, min_cells=5)

    output_dir = tmp_path / "results"
    analyzer.export(output_dir=output_dir)

    # Check that output directory was created
    assert output_dir.exists()
    assert (output_dir / "summary.txt").exists()
