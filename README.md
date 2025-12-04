# Perturbio

```
┌─────────────────────────────────────────────────────────────────┐
│                                                                 │
│  ██████╗ ███████╗██████╗ ████████╗██╗   ██╗██████╗ ██████╗ ██╗ │
│  ██╔══██╗██╔════╝██╔══██╗╚══██╔══╝██║   ██║██╔══██╗██╔══██╗██║ │
│  ██████╔╝█████╗  ██████╔╝   ██║   ██║   ██║██████╔╝██████╔╝██║ │
│  ██╔═══╝ ██╔══╝  ██╔══██╗   ██║   ██║   ██║██╔══██╗██╔══██╗██║ │
│  ██║     ███████╗██║  ██║   ██║   ╚██████╔╝██║  ██║██████╔╝██║ │
│  ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝    ╚═════╝ ╚═╝  ╚═╝╚═════╝ ╚═╝ │
│                                                                 │
│         Crop-Seq Analysis Made Simple                          │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

**Perturbio** is a comprehensive Python package for end-to-end analysis of Crop-Seq experiments (CRISPR pooled screens + single-cell RNA sequencing). From raw data to biological insights in under 5 minutes.

## Features

- **Guide Extraction**: Automatically identify CRISPR guide RNAs in single cells
- **Differential Expression**: Discover genes affected by perturbations
- **Beautiful Visualizations**: Publication-ready plots with minimal code
- **Scanpy Integration**: Works seamlessly with the scanpy ecosystem
- **Simple CLI**: One command to run complete analysis
- **Fast**: Analyze thousands of cells in minutes

## Installation

```bash
pip install perturbio
```

Or install from source:

```bash
git clone https://github.com/Perturbio/perturbio.git
cd perturbio
pip install -e .
```

## Quick Start

### Command Line

```bash
# Run complete analysis
perturbio analyze cropseq_data.h5ad --guides guides.csv

# Results saved to perturbio_results_YYYYMMDD_HHMMSS/
```

### Python API

```python
from perturbio import CropSeqAnalyzer

# One-liner magic
analyzer = CropSeqAnalyzer("cropseq_data.h5ad")
results = analyzer.run()

# Access results
print(results.top_hits("BRCA1_guide1", n=20))
results.plot_volcano("MYC_guide1")
```

### Scanpy Integration

```python
import scanpy as sc
import perturbio as pt

# Standard scanpy workflow
adata = sc.read_h5ad("cropseq_data.h5ad")
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

# Extract guides
pt.guides.extract(adata, guide_file="guides.csv")

# Differential expression
pt.tl.differential_expression(adata, groupby='perturbation', control='non-targeting')

# Visualize
pt.pl.volcano(adata, perturbation='BRCA1_guide1')
```

## Guide Library Format

Create a CSV file with your guide library:

```csv
guide_id,target_gene,guide_sequence
BRCA1_guide1,BRCA1,GCACTCAGGAAACAGCTATG
BRCA1_guide2,BRCA1,CTGAAGACTGCTCAGTGTAG
MYC_guide1,MYC,GTACTTGGTGAGGCCAGCGC
non-targeting_1,control,GTAGCGAACGTGTCCGGCGT
```

## Documentation

- [Installation Guide](docs/installation.md)
- [Quick Start Tutorial](docs/quickstart.md)
- [API Reference](https://perturbio.readthedocs.io)
- [Example Notebooks](examples/)

## Requirements

- Python 3.9+
- AnnData/Scanpy for single-cell analysis
- Works on macOS, Linux, and Windows

## Authors

This work was developed by **Siavash Ghaffari**. For any questions, feedback, or additional information, please feel free to reach out. Your input is highly valued and will help improve and refine this pipeline further.



