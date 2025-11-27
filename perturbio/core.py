"""
Core CropSeqAnalyzer class for high-level Crop-Seq analysis.
"""

from typing import Union, Optional
from pathlib import Path
import time

import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt

from perturbio.io.readers import read_h5ad, load_guide_library
from perturbio.io.writers import create_results_directory, write_summary, export_de_results, export_annotated_data
from perturbio.guides.extraction import extract
from perturbio.analysis.differential import differential_expression
from perturbio.plotting.core import volcano, umap, perturbation_counts
from perturbio.results.containers import DEResults, AnalysisResults


class CropSeqAnalyzer:
    """
    High-level interface for Crop-Seq analysis.

    This class provides a simple, user-friendly API for analyzing Crop-Seq data.
    It handles the complete workflow from loading data to generating results.

    Parameters
    ----------
    data : str, Path, or AnnData
        Path to H5AD file or AnnData object
    guide_file : str or Path, optional
        Path to guide library CSV file

    Attributes
    ----------
    adata : AnnData
        Annotated data matrix
    results : AnalysisResults
        Analysis results container

    Examples
    --------
    >>> from perturbio import CropSeqAnalyzer
    >>> # One-liner analysis
    >>> analyzer = CropSeqAnalyzer("cropseq_data.h5ad")
    >>> results = analyzer.run()
    >>>
    >>> # Step-by-step analysis
    >>> analyzer = CropSeqAnalyzer("cropseq_data.h5ad")
    >>> analyzer.extract_guides(guide_file="guides.csv")
    >>> analyzer.differential_expression()
    >>> analyzer.plot_volcano("BRCA1_guide1")
    >>> analyzer.export("results/")
    """

    def __init__(
        self,
        data: Union[str, Path, ad.AnnData],
        guide_file: Optional[Union[str, Path]] = None,
    ):
        """Initialize analyzer with data."""
        # Load data
        if isinstance(data, ad.AnnData):
            self.adata = data
            print(f"✓ Loaded AnnData: {self.adata.n_obs:,} cells × {self.adata.n_vars:,} genes")
        else:
            print("[1/5] Loading dataset...")
            self.adata = read_h5ad(data)

        # Store guide file for later
        self._guide_file = guide_file

        # Initialize results
        self.results = AnalysisResults()

        # Track what's been run
        self._guides_extracted = False
        self._de_run = False

    def extract_guides(
        self,
        guide_file: Optional[Union[str, Path]] = None,
        min_umis: int = 3,
        multiplet_threshold: float = 0.2,
    ) -> "CropSeqAnalyzer":
        """
        Extract guide barcodes from cells.

        Parameters
        ----------
        guide_file : str or Path, optional
            Path to guide library CSV. Uses value from __init__ if not provided.
        min_umis : int, default 3
            Minimum UMIs for confident guide assignment
        multiplet_threshold : float, default 0.2
            Threshold for detecting multiplets

        Returns
        -------
        CropSeqAnalyzer
            Self for method chaining
        """
        if guide_file is None:
            guide_file = self._guide_file

        if guide_file is None:
            raise ValueError(
                "Must provide guide_file either in __init__ or extract_guides()\n"
                "Example: analyzer.extract_guides(guide_file='guides.csv')"
            )

        # Extract guides
        extract(
            self.adata,
            guide_file=guide_file,
            min_umis=min_umis,
            multiplet_threshold=multiplet_threshold,
        )

        self._guides_extracted = True

        # Store guide assignments in results
        self.results.guide_assignments = self.adata.obs[[
            'guide_identity', 'perturbation', 'guide_umi_count', 'guide_confidence'
        ]].copy()

        return self

    def differential_expression(
        self,
        control_label: str = 'non-targeting',
        min_cells: int = 10,
        fdr_threshold: float = 0.05,
    ) -> "CropSeqAnalyzer":
        """
        Run differential expression analysis.

        Parameters
        ----------
        control_label : str, default 'non-targeting'
            Label for control cells
        min_cells : int, default 10
            Minimum cells required per perturbation
        fdr_threshold : float, default 0.05
            FDR threshold for significance

        Returns
        -------
        CropSeqAnalyzer
            Self for method chaining
        """
        if not self._guides_extracted:
            raise ValueError(
                "Must extract guides first\n"
                "Run: analyzer.extract_guides(guide_file='guides.csv')"
            )

        # Run DE
        differential_expression(
            self.adata,
            groupby='perturbation',
            control=control_label,
            min_cells=min_cells,
            fdr_threshold=fdr_threshold,
        )

        self._de_run = True

        # Store results
        de_df = self.adata.uns['perturbio_de']
        perturbations = de_df['perturbation'].unique().tolist()

        self.results.differential_expression = DEResults(
            data=de_df,
            perturbations=perturbations,
            control_label=control_label,
            fdr_threshold=fdr_threshold,
        )

        return self

    def run(
        self,
        guide_file: Optional[Union[str, Path]] = None,
        control_label: str = 'non-targeting',
        min_cells: int = 10,
        fdr_threshold: float = 0.05,
    ) -> AnalysisResults:
        """
        Run complete analysis pipeline.

        This is the "magic button" - runs everything:
        1. Extract guides
        2. Differential expression
        3. Return results

        Parameters
        ----------
        guide_file : str or Path, optional
            Path to guide library CSV
        control_label : str, default 'non-targeting'
            Label for control cells
        min_cells : int, default 10
            Minimum cells per perturbation
        fdr_threshold : float, default 0.05
            FDR threshold

        Returns
        -------
        AnalysisResults
            Complete analysis results

        Examples
        --------
        >>> analyzer = CropSeqAnalyzer("data.h5ad")
        >>> results = analyzer.run(guide_file="guides.csv")
        >>> print(results.summary())
        """
        start_time = time.time()

        # Extract guides
        if not self._guides_extracted:
            self.extract_guides(guide_file=guide_file)

        # Differential expression
        if not self._de_run:
            self.differential_expression(
                control_label=control_label,
                min_cells=min_cells,
                fdr_threshold=fdr_threshold,
            )

        runtime = time.time() - start_time
        print(f"\n✓ Analysis complete in {runtime//60:.0f}m {runtime%60:.0f}s")

        return self.results

    def plot_volcano(
        self,
        perturbation: str,
        fdr_threshold: float = 0.05,
        save: Optional[Union[str, Path]] = None,
    ) -> plt.Figure:
        """
        Plot volcano plot for a perturbation.

        Parameters
        ----------
        perturbation : str
            Perturbation to plot
        fdr_threshold : float, default 0.05
            FDR threshold
        save : str or Path, optional
            Path to save figure

        Returns
        -------
        Figure
            Matplotlib figure
        """
        if not self._de_run:
            raise ValueError("Must run differential_expression() first")

        return volcano(
            adata=self.adata,
            perturbation=perturbation,
            fdr_threshold=fdr_threshold,
            save=save,
        )

    def plot_umap(
        self,
        color_by: str = 'perturbation',
        save: Optional[Union[str, Path]] = None,
    ) -> plt.Figure:
        """
        Plot UMAP colored by perturbations.

        Parameters
        ----------
        color_by : str, default 'perturbation'
            Column to color by
        save : str or Path, optional
            Path to save figure

        Returns
        -------
        Figure
            Matplotlib figure
        """
        return umap(self.adata, color=color_by, save=save)

    def plot_perturbation_counts(
        self,
        save: Optional[Union[str, Path]] = None,
    ) -> plt.Figure:
        """
        Plot cell counts per perturbation.

        Parameters
        ----------
        save : str or Path, optional
            Path to save figure

        Returns
        -------
        Figure
            Matplotlib figure
        """
        if not self._guides_extracted:
            raise ValueError("Must extract guides first")

        return perturbation_counts(self.adata, save=save)

    def top_hits(self, perturbation: str, n: int = 50) -> pd.DataFrame:
        """
        Get top differentially expressed genes for a perturbation.

        Parameters
        ----------
        perturbation : str
            Perturbation label
        n : int, default 50
            Number of top genes

        Returns
        -------
        DataFrame
            Top hits
        """
        if not self._de_run:
            raise ValueError("Must run differential_expression() first")

        return self.results.top_hits(perturbation, n=n)

    def export(
        self,
        output_dir: Optional[Union[str, Path]] = None,
        formats: list = ['png'],
    ) -> Path:
        """
        Export all results to a directory.

        Parameters
        ----------
        output_dir : str or Path, optional
            Output directory. If None, creates timestamped directory.
        formats : list, default ['png']
            Figure formats to export

        Returns
        -------
        Path
            Path to results directory
        """
        if output_dir is None:
            output_dir = create_results_directory()
        else:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)

        print("\n[5/5] Exporting results...")

        # Export DE results
        if self._de_run:
            export_de_results(
                self.results.differential_expression.data,
                output_dir / "results" / "differential_expression.csv",
            )

        # Export annotated data
        if self._guides_extracted:
            export_annotated_data(
                self.adata,
                output_dir / "data" / "annotated_data.h5ad",
            )

        # Export guide assignments
        if self.results.guide_assignments is not None:
            self.results.guide_assignments.to_csv(
                output_dir / "results" / "perturbations.csv"
            )
            print(f"✓ Saved perturbations: {output_dir / 'results' / 'perturbations.csv'}")

        # Export summary
        write_summary(
            output_dir / "summary.txt",
            adata=self.adata,
            de_results=self.results.differential_expression.data if self._de_run else None,
        )

        print(f"\n✓ All results saved to: {output_dir}")
        return output_dir
