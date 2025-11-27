"""
Writers for exporting Crop-Seq analysis results.
"""

from pathlib import Path
from typing import Union, Optional, Dict, List
from datetime import datetime

import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt


def export_de_results(
    de_results: pd.DataFrame,
    output_path: Union[str, Path],
) -> None:
    """
    Export differential expression results to CSV.

    Parameters
    ----------
    de_results : DataFrame
        Differential expression results with columns:
        perturbation, gene, log_fc, pval, pval_adj, rank
    output_path : str or Path
        Path to save the CSV file

    Examples
    --------
    >>> results = analyzer.results.differential_expression
    >>> pt.io.export_de_results(results, "de_results.csv")
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    de_results.to_csv(output_path, index=False)
    print(f"✓ Saved differential expression results: {output_path}")


def export_annotated_data(
    adata: ad.AnnData,
    output_path: Union[str, Path],
) -> None:
    """
    Export AnnData object with perturbation annotations.

    Parameters
    ----------
    adata : AnnData
        Annotated AnnData object
    output_path : str or Path
        Path to save the H5AD file

    Examples
    --------
    >>> pt.io.export_annotated_data(adata, "annotated_data.h5ad")
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    adata.write_h5ad(output_path)
    print(f"✓ Saved annotated data: {output_path}")


def export_figures(
    figures: Dict[str, plt.Figure],
    output_dir: Union[str, Path],
    formats: List[str] = ['png'],
    dpi: int = 300,
) -> None:
    """
    Export matplotlib figures to files.

    Parameters
    ----------
    figures : dict
        Dictionary mapping figure names to Figure objects
    output_dir : str or Path
        Directory to save figures
    formats : list, default ['png']
        File formats to export (e.g., ['png', 'pdf', 'svg'])
    dpi : int, default 300
        Resolution for raster formats

    Examples
    --------
    >>> figures = {
    ...     'volcano_BRCA1': fig1,
    ...     'umap_perturbations': fig2,
    ... }
    >>> pt.io.export_figures(figures, "figures/", formats=['png', 'pdf'])
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for name, fig in figures.items():
        for fmt in formats:
            output_path = output_dir / f"{name}.{fmt}"
            fig.savefig(output_path, dpi=dpi, bbox_inches='tight')
            print(f"✓ Saved figure: {output_path}")


def create_results_directory(
    base_dir: Union[str, Path] = ".",
    prefix: str = "perturbio_results",
) -> Path:
    """
    Create a timestamped results directory.

    Creates a directory structure:
    perturbio_results_YYYYMMDD_HHMMSS/
    ├── data/
    ├── results/
    └── figures/
        └── volcano/

    Parameters
    ----------
    base_dir : str or Path, default "."
        Base directory for results
    prefix : str, default "perturbio_results"
        Prefix for the results directory name

    Returns
    -------
    Path
        Path to the created results directory

    Examples
    --------
    >>> output_dir = pt.io.create_results_directory()
    >>> print(output_dir)
    perturbio_results_20250126_143022
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    results_dir = Path(base_dir) / f"{prefix}_{timestamp}"

    # Create subdirectories
    (results_dir / "data").mkdir(parents=True, exist_ok=True)
    (results_dir / "results").mkdir(parents=True, exist_ok=True)
    (results_dir / "figures" / "volcano").mkdir(parents=True, exist_ok=True)

    return results_dir


def write_summary(
    output_path: Union[str, Path],
    adata: ad.AnnData,
    de_results: Optional[pd.DataFrame] = None,
    guide_stats: Optional[pd.DataFrame] = None,
    runtime: Optional[float] = None,
) -> None:
    """
    Write a human-readable summary file.

    Parameters
    ----------
    output_path : str or Path
        Path to save the summary file
    adata : AnnData
        Analyzed AnnData object
    de_results : DataFrame, optional
        Differential expression results
    guide_stats : DataFrame, optional
        Guide assignment statistics
    runtime : float, optional
        Total runtime in seconds

    Examples
    --------
    >>> pt.io.write_summary(
    ...     "summary.txt",
    ...     adata=adata,
    ...     de_results=de_results,
    ...     runtime=154.3
    ... )
    """
    output_path = Path(output_path)

    with open(output_path, 'w') as f:
        # Header
        f.write("╔" + "═" * 63 + "╗\n")
        f.write("║" + " " * 13 + "Perturbio Analysis Summary" + " " * 24 + "║\n")
        f.write("╚" + "═" * 63 + "╝\n\n")

        # Metadata
        f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"\n{'─' * 66}\n")
        f.write("DATASET OVERVIEW\n")
        f.write(f"{'─' * 66}\n")
        f.write(f"Total Cells: {adata.n_obs:,}\n")
        f.write(f"Total Genes: {adata.n_vars:,}\n\n")

        # Guide assignment
        if 'perturbation' in adata.obs.columns:
            f.write(f"{'─' * 66}\n")
            f.write("GUIDE ASSIGNMENT\n")
            f.write(f"{'─' * 66}\n")

            n_assigned = (adata.obs['perturbation'] != 'unassigned').sum()
            pct_assigned = 100 * n_assigned / adata.n_obs

            f.write(f"Cells with guides: {n_assigned:,} ({pct_assigned:.1f}%)\n")

            # Count controls
            perturbations = adata.obs['perturbation'].value_counts()
            n_control = perturbations.get('non-targeting', 0)
            if n_control > 0:
                pct_control = 100 * n_control / adata.n_obs
                f.write(f"Non-targeting controls: {n_control:,} ({pct_control:.1f}%)\n")

            n_unassigned = (adata.obs['perturbation'] == 'unassigned').sum()
            if n_unassigned > 0:
                pct_unassigned = 100 * n_unassigned / adata.n_obs
                f.write(f"Unassigned cells: {n_unassigned:,} ({pct_unassigned:.1f}%)\n")

            f.write(f"\nGuides tested: {len(perturbations):,}\n")

            if len(perturbations) > 0:
                median_cells = perturbations.median()
                f.write(f"Cells per guide (median): {median_cells:.0f}\n")

        # Differential expression
        if de_results is not None:
            f.write(f"\n{'─' * 66}\n")
            f.write("DIFFERENTIAL EXPRESSION\n")
            f.write(f"{'─' * 66}\n")

            n_sig = (de_results['pval_adj'] < 0.05).sum()
            f.write(f"Total significant genes: {n_sig:,} (FDR < 0.05)\n")

            # Top perturbations
            if 'perturbation' in de_results.columns:
                top_perturbs = de_results.groupby('perturbation').apply(
                    lambda x: (x['pval_adj'] < 0.05).sum()
                ).nlargest(5)

                f.write(f"\nTop 5 Perturbations by Number of Significant Genes:\n")
                for i, (pert, n_genes) in enumerate(top_perturbs.items(), 1):
                    # Get top gene for this perturbation
                    pert_data = de_results[de_results['perturbation'] == pert].iloc[0]
                    gene = pert_data['gene']
                    log_fc = pert_data['log_fc']
                    pval_adj = pert_data['pval_adj']

                    f.write(
                        f"{i}. {pert:<15} - {n_genes:>3} sig genes  - "
                        f"{gene} ↓{abs(log_fc):.1f} (FDR={pval_adj:.1e})\n"
                    )

        # Runtime
        if runtime is not None:
            f.write(f"\n{'─' * 66}\n")
            f.write("OUTPUT FILES\n")
            f.write(f"{'─' * 66}\n")
            f.write("✓ Differential expression results\n")
            f.write("✓ Annotated dataset with perturbations\n")
            f.write("✓ Visualization plots\n")
            f.write("✓ Summary statistics\n")
            f.write(f"\nTotal analysis time: {runtime//60:.0f}m {runtime%60:.0f}s\n")

    print(f"✓ Saved summary: {output_path}")
