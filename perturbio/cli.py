"""
Command-line interface for Perturbio.
"""

import sys
from pathlib import Path

import click

from perturbio import __version__
from perturbio.core import CropSeqAnalyzer


@click.group()
@click.version_option(version=__version__)
def main():
    """
    Perturbio - Crop-Seq Analysis Made Simple

    Analyze CRISPR pooled screens with single-cell RNA-seq.
    """
    pass


@main.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.option(
    '--guides',
    type=click.Path(exists=True),
    help='Path to guide library CSV file'
)
@click.option(
    '--output',
    '-o',
    type=click.Path(),
    default=None,
    help='Output directory (default: perturbio_results_TIMESTAMP)'
)
@click.option(
    '--control',
    default='non-targeting',
    help='Control cell label (default: non-targeting)'
)
@click.option(
    '--min-cells',
    default=10,
    type=int,
    help='Minimum cells per perturbation (default: 10)'
)
@click.option(
    '--min-umis',
    default=3,
    type=int,
    help='Minimum UMIs for guide assignment (default: 3)'
)
@click.option(
    '--fdr',
    default=0.05,
    type=float,
    help='FDR threshold for significance (default: 0.05)'
)
def analyze(input_file, guides, output, control, min_cells, min_umis, fdr):
    """
    Run complete Crop-Seq analysis pipeline.

    INPUT_FILE: Path to H5AD file with single-cell data

    Examples:

        # Basic usage
        perturbio analyze cropseq_data.h5ad --guides guides.csv

        # With custom parameters
        perturbio analyze data.h5ad --guides guides.csv --output results/ --min-cells 20
    """
    click.echo("┌" + "─" * 63 + "┐")
    click.echo("│ " + click.style("Perturbio", fg="blue", bold=True) + " v" + __version__ + " - Crop-Seq Analysis Pipeline" + " " * 18 + "│")
    click.echo("└" + "─" * 63 + "┘\n")

    try:
        # Initialize analyzer
        analyzer = CropSeqAnalyzer(input_file, guide_file=guides)

        # Run analysis
        analyzer.run(
            guide_file=guides,
            control_label=control,
            min_cells=min_cells,
            fdr_threshold=fdr,
        )

        # Generate visualizations
        click.echo("\n[4/5] Generating visualizations...")

        # Try to create plots (may fail if UMAP not computed)
        try:
            analyzer.plot_perturbation_counts()
            click.echo("  ✓ Perturbation assignment counts")
        except Exception as e:
            click.echo(f"  ⚠ Skipping perturbation counts: {str(e)}", err=True)

        try:
            # Check if UMAP exists
            if 'X_umap' in analyzer.adata.obsm:
                analyzer.plot_umap(color_by='perturbation')
                click.echo("  ✓ UMAP colored by perturbation")
            else:
                click.echo("  ⚠ Skipping UMAP (run sc.tl.umap first)")
        except Exception as e:
            click.echo(f"  ⚠ Skipping UMAP: {str(e)}", err=True)

        # Create volcano plots for top perturbations
        try:
            top_perts = analyzer.results.differential_expression.perturbations[:5]
            for pert in top_perts:
                analyzer.plot_volcano(pert)
            click.echo(f"  ✓ Volcano plots ({len(top_perts)} perturbations)")
        except Exception as e:
            click.echo(f"  ⚠ Skipping volcano plots: {str(e)}", err=True)

        # Export results
        output_dir = analyzer.export(output_dir=output)

        # Final summary
        click.echo("\n┌" + "─" * 63 + "┐")
        click.echo("│" + " " * 20 + click.style("✓ Analysis Complete!", fg="green", bold=True) + " " * 20 + "│")
        click.echo("├" + "─" * 63 + "┤")
        click.echo(f"│ Results saved to: {str(output_dir):<44} │")
        click.echo("│" + " " * 63 + "│")
        click.echo("│ Next steps:" + " " * 52 + "│")
        click.echo("│  • Review summary: results/summary.txt" + " " * 24 + "│")
        click.echo("│  • Check top hits: results/differential_expression.csv" + " " * 8 + "│")
        click.echo("│  • View volcano plots: figures/volcano/*.png" + " " * 18 + "│")
        click.echo("└" + "─" * 63 + "┘")

    except Exception as e:
        click.echo("\n┌" + "─" * 63 + "┐", err=True)
        click.echo("│ " + click.style("✗ Error", fg="red", bold=True) + " " * 56 + "│", err=True)
        click.echo("├" + "─" * 63 + "┤", err=True)
        click.echo(f"│ {str(e):<61} │", err=True)
        click.echo("└" + "─" * 63 + "┘", err=True)
        sys.exit(1)


@main.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.option('--guides', type=click.Path(exists=True), required=True, help='Guide library CSV')
@click.option('--output', '-o', type=click.Path(), required=True, help='Output H5AD file')
@click.option('--min-umis', default=3, type=int, help='Minimum UMIs (default: 3)')
def extract_guides(input_file, guides, output, min_umis):
    """
    Extract CRISPR guide barcodes from cells.

    INPUT_FILE: Path to H5AD file with single-cell data
    """
    try:
        analyzer = CropSeqAnalyzer(input_file)
        analyzer.extract_guides(guide_file=guides, min_umis=min_umis)

        # Save annotated data
        analyzer.adata.write_h5ad(output)
        click.echo(f"\n✓ Saved annotated data: {output}")

    except Exception as e:
        click.echo(f"\n✗ Error: {str(e)}", err=True)
        sys.exit(1)


@main.command()
def version():
    """Show Perturbio version."""
    click.echo(f"Perturbio version {__version__}")


if __name__ == '__main__':
    main()
