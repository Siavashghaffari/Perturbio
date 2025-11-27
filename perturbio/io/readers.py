"""
Readers for loading Crop-Seq data files.
"""

from pathlib import Path
from typing import Union, Optional
import warnings

import anndata as ad
import pandas as pd


def read_h5ad(
    file_path: Union[str, Path],
    validate: bool = True,
) -> ad.AnnData:
    """
    Read an AnnData object from an H5AD file.

    Parameters
    ----------
    file_path : str or Path
        Path to the H5AD file
    validate : bool, default True
        Whether to validate the loaded data

    Returns
    -------
    AnnData
        Loaded AnnData object

    Raises
    ------
    FileNotFoundError
        If the file does not exist
    ValueError
        If validation fails and the file is not a valid AnnData object

    Examples
    --------
    >>> import perturbio as pt
    >>> adata = pt.io.read_h5ad("cropseq_data.h5ad")
    >>> print(adata)
    """
    file_path = Path(file_path)

    # Check file exists
    if not file_path.exists():
        raise FileNotFoundError(
            f"Could not find file: {file_path}\n"
            f"Suggestions:\n"
            f"  • Check that the file path is correct\n"
            f"  • Ensure the file has .h5ad extension\n"
            f"  • Try using an absolute path"
        )

    # Check file extension
    if file_path.suffix != ".h5ad":
        warnings.warn(
            f"File does not have .h5ad extension: {file_path.suffix}\n"
            f"Attempting to read anyway...",
            UserWarning
        )

    # Load the data
    try:
        adata = ad.read_h5ad(file_path)
    except Exception as e:
        raise ValueError(
            f"Failed to read H5AD file: {file_path}\n"
            f"Error: {str(e)}\n"
            f"Ensure the file is a valid AnnData H5AD file"
        ) from e

    # Validate if requested
    if validate:
        validate_adata(adata)

    return adata


def validate_adata(adata: ad.AnnData) -> None:
    """
    Validate an AnnData object for Crop-Seq analysis.

    Checks:
    - Has cell × gene matrix
    - Contains expression data
    - Has proper dimensions

    Parameters
    ----------
    adata : AnnData
        AnnData object to validate

    Raises
    ------
    ValueError
        If validation fails

    Examples
    --------
    >>> import perturbio as pt
    >>> import scanpy as sc
    >>> adata = sc.read_h5ad("data.h5ad")
    >>> pt.io.validate_adata(adata)
    """
    # Check that we have data
    if adata.n_obs == 0:
        raise ValueError(
            "AnnData object has no cells (n_obs = 0)\n"
            "Cannot analyze empty dataset"
        )

    if adata.n_vars == 0:
        raise ValueError(
            "AnnData object has no genes (n_vars = 0)\n"
            "Cannot analyze empty dataset"
        )

    # Check for expression matrix
    if adata.X is None:
        # Look for data in layers
        if len(adata.layers) > 0:
            layer_names = ", ".join(adata.layers.keys())
            raise ValueError(
                "Expected count matrix in adata.X, but found None\n"
                f"Your file contains layers: [{layer_names}]\n"
                f"Solutions:\n"
                f"  • Copy a layer to X: adata.X = adata.layers['raw']\n"
                f"  • Specify which layer to use when loading"
            )
        else:
            raise ValueError(
                "Expected count matrix in adata.X, but found None\n"
                "AnnData object has no expression data"
            )

    # Warn about data normalization
    if adata.X.max() > 100:
        warnings.warn(
            f"Data appears unnormalized (max value: {adata.X.max():,.0f})\n"
            "Consider normalizing with scanpy:\n"
            "  sc.pp.normalize_total(adata)\n"
            "  sc.pp.log1p(adata)",
            UserWarning
        )

    # Info about dataset
    print(f"✓ Loaded valid AnnData object")
    print(f"  • {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    if adata.X.max() <= 20:
        print(f"  • Detected normalized/log-transformed counts")
    else:
        print(f"  • Raw counts (max value: {adata.X.max():,.0f})")


def load_guide_library(
    file_path: Union[str, Path],
    required_columns: Optional[list] = None,
) -> pd.DataFrame:
    """
    Load a guide library from a CSV file.

    Expected format:
    - guide_id: Unique identifier for each guide
    - target_gene: Gene being targeted (or "control")
    - guide_sequence: Optional RNA sequence

    Parameters
    ----------
    file_path : str or Path
        Path to the guide library CSV file
    required_columns : list, optional
        List of required column names. If None, uses ['guide_id', 'target_gene']

    Returns
    -------
    DataFrame
        Guide library with validated columns

    Raises
    ------
    FileNotFoundError
        If file does not exist
    ValueError
        If required columns are missing

    Examples
    --------
    >>> import perturbio as pt
    >>> guides = pt.io.load_guide_library("guides.csv")
    >>> print(guides.head())
    """
    if required_columns is None:
        required_columns = ['guide_id', 'target_gene']

    file_path = Path(file_path)

    # Check file exists
    if not file_path.exists():
        raise FileNotFoundError(
            f"Could not find guide library: {file_path}\n"
            f"Ensure the file path is correct"
        )

    # Load CSV
    try:
        guides = pd.read_csv(file_path)
    except Exception as e:
        raise ValueError(
            f"Failed to read CSV file: {file_path}\n"
            f"Error: {str(e)}"
        ) from e

    # Validate required columns
    missing_cols = [col for col in required_columns if col not in guides.columns]
    if missing_cols:
        raise ValueError(
            f"Guide library missing required columns: {missing_cols}\n"
            f"Found columns: {list(guides.columns)}\n"
            f"Required columns: {required_columns}\n\n"
            f"Expected format:\n"
            f"  guide_id,target_gene,guide_sequence\n"
            f"  BRCA1_guide1,BRCA1,GCACTCAGGAAACAGCTATG\n"
            f"  MYC_guide1,MYC,GTACTTGGTGAGGCCAGCGC"
        )

    # Check for duplicates
    if guides['guide_id'].duplicated().any():
        duplicates = guides[guides['guide_id'].duplicated()]['guide_id'].tolist()
        warnings.warn(
            f"Found duplicate guide_ids: {duplicates}\n"
            "This may cause issues with guide assignment",
            UserWarning
        )

    print(f"✓ Loaded guide library: {len(guides)} guides")

    # Count targeting vs control guides
    n_control = guides['target_gene'].str.contains('control|non-targeting|NTC', case=False, na=False).sum()
    n_targeting = len(guides) - n_control
    print(f"  • {n_targeting} targeting guides")
    print(f"  • {n_control} control guides")

    return guides
