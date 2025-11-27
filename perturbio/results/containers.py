"""
Container classes for storing and accessing analysis results.
"""

from dataclasses import dataclass
from typing import Optional, List
import pandas as pd


@dataclass
class DEResults:
    """
    Container for differential expression results.

    Attributes
    ----------
    data : DataFrame
        Full differential expression results with columns:
        perturbation, gene, log_fc, pval, pval_adj, rank
    perturbations : list
        List of tested perturbations
    control_label : str
        Label used for control cells
    fdr_threshold : float
        FDR threshold used for significance

    Examples
    --------
    >>> results = DEResults(de_df, perturbations=['BRCA1', 'MYC'], control_label='non-targeting')
    >>> top = results.top_hits('BRCA1', n=20)
    >>> sig = results.significant_genes(fdr_threshold=0.05)
    """
    data: pd.DataFrame
    perturbations: List[str]
    control_label: str
    fdr_threshold: float = 0.05

    def top_hits(self, perturbation: str, n: int = 50) -> pd.DataFrame:
        """
        Get top differentially expressed genes for a perturbation.

        Parameters
        ----------
        perturbation : str
            Perturbation label
        n : int, default 50
            Number of top genes to return

        Returns
        -------
        DataFrame
            Top n genes ranked by significance
        """
        if perturbation not in self.perturbations:
            raise ValueError(
                f"Perturbation '{perturbation}' not found\n"
                f"Available: {self.perturbations}"
            )

        pert_data = self.data[self.data['perturbation'] == perturbation]
        return pert_data.head(n)[['gene', 'log_fc', 'pval', 'pval_adj']]

    def significant_genes(self, fdr_threshold: Optional[float] = None) -> pd.DataFrame:
        """
        Get all significantly perturbed genes across all perturbations.

        Parameters
        ----------
        fdr_threshold : float, optional
            FDR threshold. If None, uses self.fdr_threshold

        Returns
        -------
        DataFrame
            Significant genes
        """
        if fdr_threshold is None:
            fdr_threshold = self.fdr_threshold

        sig = self.data[self.data['pval_adj'] < fdr_threshold]
        return sig.sort_values('pval_adj')

    def genes_per_perturbation(self, fdr_threshold: Optional[float] = None) -> pd.Series:
        """
        Count significant genes for each perturbation.

        Parameters
        ----------
        fdr_threshold : float, optional
            FDR threshold

        Returns
        -------
        Series
            Number of significant genes per perturbation
        """
        if fdr_threshold is None:
            fdr_threshold = self.fdr_threshold

        sig = self.data[self.data['pval_adj'] < fdr_threshold]
        return sig.groupby('perturbation').size().sort_values(ascending=False)


@dataclass
class AnalysisResults:
    """
    Container for complete Crop-Seq analysis results.

    Attributes
    ----------
    differential_expression : DEResults
        Differential expression results
    guide_assignments : DataFrame
        Guide assignments per cell
    perturbations_summary : DataFrame
        Summary statistics per perturbation

    Examples
    --------
    >>> results = analyzer.results
    >>> print(results.differential_expression.top_hits('BRCA1'))
    >>> print(results.perturbations_summary)
    """
    differential_expression: Optional[DEResults] = None
    guide_assignments: Optional[pd.DataFrame] = None
    perturbations_summary: Optional[pd.DataFrame] = None

    def summary(self) -> str:
        """
        Get a text summary of the analysis results.

        Returns
        -------
        str
            Formatted summary
        """
        lines = ["Crop-Seq Analysis Results", "=" * 40]

        if self.guide_assignments is not None:
            n_cells = len(self.guide_assignments)
            n_assigned = (self.guide_assignments['perturbation'] != 'unassigned').sum()
            lines.append(f"\nCells: {n_cells:,}")
            lines.append(f"Assigned: {n_assigned:,} ({100*n_assigned/n_cells:.1f}%)")

        if self.differential_expression is not None:
            n_pert = len(self.differential_expression.perturbations)
            n_sig = len(self.differential_expression.significant_genes())
            lines.append(f"\nPerturbations tested: {n_pert}")
            lines.append(f"Significant genes: {n_sig:,}")

        return "\n".join(lines)

    def top_hits(self, perturbation: str, n: int = 50) -> pd.DataFrame:
        """
        Convenience method to get top hits for a perturbation.

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
        if self.differential_expression is None:
            raise ValueError("No differential expression results available")

        return self.differential_expression.top_hits(perturbation, n=n)
