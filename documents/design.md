# Perturbio - Interface Design

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

---

## Table of Contents

1. [CLI Command Structure](#cli-command-structure)
2. [Terminal Output Examples](#terminal-output-examples)
3. [Python API Design](#python-api-design)
4. [Package Architecture](#package-architecture)
5. [User Journey Flows](#user-journey-flows)
6. [File Structure & Outputs](#file-structure--outputs)

---

## CLI Command Structure

### Command Hierarchy

```
perturbio
│
├── analyze              # Main analysis command (MVP)
│   ├── --input         # Path to h5ad or BAM file
│   ├── --guides        # Path to guide library CSV
│   ├── --output        # Output directory
│   ├── --control       # Control label (default: auto-detect)
│   ├── --min-cells     # Minimum cells per perturbation
│   └── --fdr           # FDR threshold for significance
│
├── extract-guides       # Standalone guide extraction
│   ├── --input
│   ├── --guides
│   └── --output
│
├── diff-exp            # Standalone differential expression
│   ├── --input
│   ├── --groupby
│   └── --control
│
├── plot                # Generate plots from analyzed data
│   ├── --input
│   ├── --type          # volcano, umap, counts
│   └── --perturbation
│
├── version             # Show version
└── help                # Show help
```

### Usage Examples

```bash
# Quickstart (all defaults)
$ perturbio analyze cropseq_data.h5ad

# With custom guide library
$ perturbio analyze \
    --input cropseq_data.h5ad \
    --guides my_guides.csv \
    --output results/

# Advanced options
$ perturbio analyze \
    --input cropseq_data.h5ad \
    --guides guides.csv \
    --output results/ \
    --control "non-targeting" \
    --min-cells 10 \
    --fdr 0.05

# Standalone guide extraction
$ perturbio extract-guides \
    --input cropseq_data.h5ad \
    --guides guides.csv \
    --output annotated_data.h5ad

# Generate specific plot
$ perturbio plot \
    --input analyzed_data.h5ad \
    --type volcano \
    --perturbation BRCA1_guide1
```

---

## Terminal Output Examples

### Example 1: Successful Analysis Run

```
$ perturbio analyze cropseq_data.h5ad --guides guides.csv

┌─────────────────────────────────────────────────────────────────┐
│ Perturbio v0.1.0 - Crop-Seq Analysis Pipeline                  │
└─────────────────────────────────────────────────────────────────┘

[1/5] Loading dataset...
  ✓ Loaded cropseq_data.h5ad
  ✓ 8,467 cells × 23,451 genes
  ✓ Detected normalized counts in layer 'X'

[2/5] Extracting CRISPR guide barcodes...
  ✓ Loaded guide library: 25 guides (24 targeting + 1 control)
  ↻ Assigning guides to cells... ━━━━━━━━━━━━━━━━━━━━━━ 100%

  Guide Assignment Summary:
  ┌──────────────────────────┬────────┐
  │ Category                 │ Count  │
  ├──────────────────────────┼────────┤
  │ Cells with guides        │  8,102 │
  │ Non-targeting controls   │    289 │
  │ Unassigned cells         │    76  │
  │ Low confidence           │     0  │
  └──────────────────────────┴────────┘

[3/5] Running differential expression...
  ↻ Testing 24 perturbations vs control... ━━━━━━━━━━━━━ 100%

  Perturbation Results:
  ┌─────────────────┬───────┬──────────┬─────────────┐
  │ Perturbation    │ Cells │ Sig Gene │ Top Gene    │
  ├─────────────────┼───────┼──────────┼─────────────┤
  │ BRCA1_guide1    │  347  │   421    │ BRCA1 ↓3.8  │
  │ MYC_guide1      │  298  │   567    │ MYC ↓4.2    │
  │ TP53_guide1     │  412  │   689    │ TP53 ↓3.1   │
  │ EGFR_guide2     │  221  │   234    │ EGFR ↓2.9   │
  │ ...             │  ...  │   ...    │ ...         │
  └─────────────────┴───────┴──────────┴─────────────┘

  ✓ Found 3,847 significantly perturbed genes (FDR < 0.05)

[4/5] Generating visualizations...
  ✓ Perturbation assignment counts
  ✓ UMAP colored by perturbation
  ✓ Volcano plots (24 perturbations)

[5/5] Exporting results...
  ✓ Saved differential expression results
  ✓ Saved annotated data with perturbations
  ✓ Saved figures (27 plots)

┌─────────────────────────────────────────────────────────────────┐
│                    ✓ Analysis Complete!                         │
├─────────────────────────────────────────────────────────────────┤
│ Results saved to: ./perturbio_results_20250126_143022/         │
│                                                                 │
│ Next steps:                                                     │
│  • Review summary: results/summary.txt                          │
│  • Check top hits: results/top_hits_summary.csv                 │
│  • View volcano plots: results/figures/volcano_*.png            │
│                                                                 │
│ Total time: 2m 34s                                              │
└─────────────────────────────────────────────────────────────────┘
```

### Example 2: Analysis with Warnings

```
$ perturbio analyze cropseq_data.h5ad

[1/5] Loading dataset...
  ✓ Loaded cropseq_data.h5ad
  ⚠ Warning: Data appears unnormalized (max value: 45,821)
  → Applying log-normalization automatically

[2/5] Extracting CRISPR guide barcodes...
  ⚠ Warning: No guide library provided
  → Auto-detecting guides from gene names
  ✓ Found 18 potential guide sequences

  Guide Assignment Summary:
  ┌──────────────────────────┬────────┐
  │ Category                 │ Count  │
  ├──────────────────────────┼────────┤
  │ Cells with guides        │  5,234 │
  │ Non-targeting controls   │    142 │
  │ Unassigned cells         │  1,876 │
  │ Low confidence           │     48 │
  └──────────────────────────┴────────┘

  ⚠ Warning: 1,876 cells (31%) have no guide detected
  → Consider adjusting --min-umis threshold

[3/5] Running differential expression...
  ⚠ Warning: KRAS_guide2 has only 8 cells (min: 10)
  → Skipping KRAS_guide2 from analysis
  ↻ Testing 17 perturbations vs control... ━━━━━━━━━━━━━ 100%

...
```

### Example 3: Error Handling

```
$ perturbio analyze missing_file.h5ad

┌─────────────────────────────────────────────────────────────────┐
│ ✗ Error: File not found                                        │
├─────────────────────────────────────────────────────────────────┤
│ Could not find: missing_file.h5ad                               │
│                                                                 │
│ Suggestions:                                                    │
│  • Check that the file path is correct                          │
│  • Ensure the file has .h5ad extension                          │
│  • Try using an absolute path                                   │
│                                                                 │
│ For help: perturbio --help                                      │
└─────────────────────────────────────────────────────────────────┘
```

```
$ perturbio analyze cropseq_data.h5ad --guides guides.csv

[1/5] Loading dataset...
  ✓ Loaded cropseq_data.h5ad

[2/5] Extracting CRISPR guide barcodes...
  ✓ Loaded guide library: 25 guides
  ✗ Error: No guides detected in dataset

┌─────────────────────────────────────────────────────────────────┐
│ ✗ Guide Extraction Failed                                      │
├─────────────────────────────────────────────────────────────────┤
│ No genes in the dataset matched the guide library.              │
│                                                                 │
│ Your guide library contains:                                    │
│  • BRCA1_guide1, BRCA1_guide2, MYC_guide1, ...                  │
│                                                                 │
│ Your dataset contains genes like:                               │
│  • ENSG00000139618, ENSG00000136997, ...                        │
│                                                                 │
│ Possible issues:                                                │
│  ✗ Guide names don't match gene IDs in dataset                  │
│  ✗ Guides may be in a separate assay (check .obsm or .layers)  │
│  ✗ Wrong guide library file                                     │
│                                                                 │
│ Solutions:                                                      │
│  • Ensure guide names in CSV match gene names in h5ad           │
│  • Convert gene IDs to symbols, or vice versa                   │
│  • Check data preprocessing steps                               │
│                                                                 │
│ For help: https://perturbio.readthedocs.io/troubleshooting     │
└─────────────────────────────────────────────────────────────────┘
```

---

## Python API Design

### High-Level API (Beginner-Friendly)

```python
"""
The "Magic Button" API - everything in 2 lines
"""
from perturbio import CropSeqAnalyzer

# One-liner analysis
analyzer = CropSeqAnalyzer("cropseq_data.h5ad")
results = analyzer.run()  # Extract guides → DE → plots → export

# Access results
print(results.summary())
print(results.top_hits("BRCA1_guide1", n=20))
results.plot_volcano("MYC_guide1")
```

### Step-by-Step API (More Control)

```python
"""
Step-by-step workflow for custom analysis
"""
from perturbio import CropSeqAnalyzer

# Initialize
analyzer = CropSeqAnalyzer("cropseq_data.h5ad")

# Step 1: Extract guides
analyzer.extract_guides(
    guide_file="guides.csv",
    min_umis=3,
    control_label="non-targeting"
)

# Check guide assignment
print(analyzer.guide_summary())
# Output:
#   Cells with guides: 8,102 (95.7%)
#   Non-targeting controls: 289 (3.4%)
#   Unassigned: 76 (0.9%)

# Step 2: Differential expression
analyzer.differential_expression(
    control_label="non-targeting",
    fdr_threshold=0.05,
    min_cells=10
)

# View results
de_results = analyzer.results.differential_expression
print(de_results.head())

# Step 3: Visualize
analyzer.plot_perturbation_counts()
analyzer.plot_umap(color_by="perturbation")
analyzer.plot_volcano("BRCA1_guide1", top_n=20)

# Step 4: Export
analyzer.export(output_dir="results/")
```

### Low-Level API (Scanpy Integration)

```python
"""
Direct integration with scanpy workflows
"""
import scanpy as sc
import perturbio as pt

# Standard scanpy preprocessing
adata = sc.read_h5ad("cropseq_data.h5ad")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Extract guides (modifies adata in-place)
pt.guides.extract(
    adata,
    guide_file="guides.csv",
    key_added="perturbation"
)

# Check what was added
print(adata.obs.columns)
# [..., 'perturbation', 'guide_identity', 'guide_umi_count']

# Dimensionality reduction
sc.pp.highly_variable_genes(adata)
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

# Visualize perturbations
sc.pl.umap(adata, color='perturbation', legend_loc='on data')

# Differential expression the Perturbio way
pt.tl.differential_expression(
    adata,
    groupby='perturbation',
    control='non-targeting',
    key_added='perturbio_de'
)

# Results stored in adata.uns
de_results = adata.uns['perturbio_de']

# Generate volcano plot
pt.pl.volcano(
    adata,
    perturbation='BRCA1_guide1',
    de_key='perturbio_de'
)
```

### Results Object API

```python
"""
Working with analysis results
"""
# Access results
results = analyzer.results

# Differential expression results
results.differential_expression
# Returns: pandas DataFrame with columns:
#   perturbation | gene | log_fc | pval | pval_adj | rank

# Get top hits for a perturbation
top_genes = results.top_hits("BRCA1_guide1", n=50)
print(top_genes)
#   gene      log_fc  pval_adj
#   BRCA1     -3.82   1.2e-45
#   RAD51     -2.14   3.4e-23
#   BRIP1     -1.89   7.8e-18

# Get all significant genes across all perturbations
sig_genes = results.significant_genes(fdr_threshold=0.05)
print(f"Found {len(sig_genes)} significant genes")

# Perturbation summary
print(results.perturbations_summary)
#   perturbation    cells  sig_genes  top_gene
#   BRCA1_guide1     347      421      BRCA1
#   MYC_guide1       298      567      MYC

# Guide assignment info
print(results.guide_assignments)
#   cell_barcode      guide_identity  perturbation  umi_count
#   AAACCTGAGCGATGAC  BRCA1_guide1    BRCA1         12
#   AAACCTGCACGGTGTC  non-targeting   control       5

# Export specific results
results.export_top_hits("top_hits.csv", top_n=100)
results.export_de_results("all_de_results.csv")
```

---

## Package Architecture

### High-Level Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                         User Interface                          │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  ┌──────────────────┐              ┌──────────────────┐        │
│  │   CLI (Click)    │              │   Python API     │        │
│  │                  │              │                  │        │
│  │  $ perturbio     │              │  from perturbio  │        │
│  │    analyze       │              │  import *        │        │
│  └────────┬─────────┘              └─────────┬────────┘        │
│           │                                  │                 │
│           └──────────────┬───────────────────┘                 │
│                          │                                     │
└──────────────────────────┼─────────────────────────────────────┘
                           │
┌──────────────────────────┼─────────────────────────────────────┐
│                          ▼                                     │
│                  CropSeqAnalyzer                                │
│              (High-Level Orchestrator)                          │
│                                                                 │
│  Methods:                                                       │
│   • run() - Full pipeline                                       │
│   • extract_guides()                                            │
│   • differential_expression()                                   │
│   • plot_*()                                                    │
│   • export()                                                    │
└──────────────────────────┬─────────────────────────────────────┘
                           │
        ┌──────────────────┼──────────────────┐
        │                  │                  │
        ▼                  ▼                  ▼
┌───────────────┐  ┌──────────────┐  ┌──────────────┐
│  io.readers   │  │    guides    │  │   analysis   │
│               │  │   .extract   │  │ .differential│
│ • read_h5ad() │  │              │  │              │
│ • read_bam()  │  │ • assign()   │  │ • wilcoxon() │
│ • validate()  │  │ • detect()   │  │ • rank_genes │
└───────────────┘  └──────────────┘  └──────────────┘
        │                  │                  │
        └──────────────────┼──────────────────┘
                           │
                           ▼
                  ┌─────────────────┐
                  │    AnnData      │
                  │   (scanpy)      │
                  │                 │
                  │  .obs           │
                  │  .var           │
                  │  .X             │
                  │  .uns           │
                  └─────────────────┘
                           │
        ┌──────────────────┼──────────────────┐
        │                  │                  │
        ▼                  ▼                  ▼
┌───────────────┐  ┌──────────────┐  ┌──────────────┐
│   plotting    │  │    results   │  │   io.export  │
│               │  │              │  │              │
│ • volcano()   │  │ • DEResults  │  │ • to_csv()   │
│ • umap()      │  │ • top_hits() │  │ • to_h5ad()  │
│ • counts()    │  │ • summary()  │  │ • figures()  │
└───────────────┘  └──────────────┘  └──────────────┘
```

### Module Structure

```
perturbio/
│
├── __init__.py                  # Public API exports
│   └── Exports: CropSeqAnalyzer, extract_guides, plot_volcano, etc.
│
├── core.py                      # Main analyzer class
│   └── class CropSeqAnalyzer
│       ├── __init__(data: str | AnnData)
│       ├── extract_guides() -> Self
│       ├── differential_expression() -> Self
│       ├── run() -> Self
│       ├── plot_volcano()
│       ├── plot_umap()
│       ├── plot_perturbation_counts()
│       └── export()
│
├── io/
│   ├── __init__.py
│   ├── readers.py               # Load data files
│   │   ├── read_h5ad()
│   │   ├── read_10x_h5()
│   │   └── validate_adata()
│   └── writers.py               # Export results
│       ├── export_de_results()
│       ├── export_annotated_data()
│       └── export_figures()
│
├── guides/
│   ├── __init__.py
│   ├── extraction.py            # Guide detection
│   │   ├── extract_guides()
│   │   ├── assign_guides_to_cells()
│   │   ├── detect_multiplets()
│   │   └── filter_low_confidence()
│   └── library.py               # Guide library handling
│       ├── load_guide_library()
│       └── validate_guide_library()
│
├── analysis/
│   ├── __init__.py
│   └── differential.py          # DE analysis
│       ├── differential_expression()
│       ├── wilcoxon_test()
│       └── rank_genes_per_perturbation()
│
├── plotting/
│   ├── __init__.py
│   └── core.py                  # Visualization functions
│       ├── volcano()
│       ├── umap()
│       ├── perturbation_counts()
│       └── _style_plot()        # Internal styling
│
├── results/
│   ├── __init__.py
│   └── containers.py            # Results data structures
│       ├── class DEResults
│       ├── class GuideAssignments
│       └── class AnalysisResults
│
├── utils/
│   ├── __init__.py
│   ├── logging.py               # Logging setup
│   └── validation.py            # Input validation
│
└── cli.py                       # Command-line interface
    ├── main()
    ├── analyze()
    ├── extract_guides()
    └── plot()
```

### Data Flow Diagram

```
┌─────────────┐
│   Input     │
│  .h5ad      │
│   File      │
└──────┬──────┘
       │
       ▼
┌─────────────────────────────────────┐
│  1. Load & Validate                 │
│     io.readers.read_h5ad()          │
│     ✓ Check for count matrix        │
│     ✓ Validate dimensions           │
└──────┬──────────────────────────────┘
       │
       ▼  AnnData(8467 cells × 23451 genes)
       │
┌─────────────────────────────────────┐
│  2. Extract Guide Barcodes          │
│     guides.extract_guides()         │
│     • Match guides to genes         │
│     • Assign cells to guides        │
│     • Flag multiplets/low quality   │
└──────┬──────────────────────────────┘
       │
       ▼  AnnData + adata.obs['perturbation']
       │
┌─────────────────────────────────────┐
│  3. Differential Expression         │
│     analysis.differential()         │
│     • For each perturbation:        │
│       - Compare to control          │
│       - Wilcoxon test               │
│       - Compute log FC              │
│       - Adjust p-values (FDR)       │
└──────┬──────────────────────────────┘
       │
       ▼  DEResults(24 perturbations × 23451 genes)
       │
┌─────────────────────────────────────┐
│  4. Visualization                   │
│     plotting.volcano()              │
│     plotting.umap()                 │
│     plotting.perturbation_counts()  │
└──────┬──────────────────────────────┘
       │
       ▼
┌─────────────────────────────────────┐
│  5. Export Results                  │
│     io.writers.export_*()           │
│     • CSV files (DE results)        │
│     • Annotated h5ad                │
│     • Figure files (PNG/PDF)        │
└──────┬──────────────────────────────┘
       │
       ▼
┌─────────────┐
│  Output     │
│  Directory  │
│  results/   │
└─────────────┘
```

---

## User Journey Flows

### Journey 1: First-Time User (Quickstart)

```
┌─────────────────────────────────────────────────────────────────┐
│ Persona: Graduate student, familiar with scanpy                 │
│ Goal: Analyze Crop-Seq data for the first time                  │
│ Time: 5 minutes                                                  │
└─────────────────────────────────────────────────────────────────┘

Step 1: Install
┌────────────────────────────────────┐
│ $ pip install perturbio            │
│                                    │
│ Installing...                      │
│ ✓ Successfully installed perturbio │
└────────────────────────────────────┘
           │
           ▼
Step 2: Prepare data
┌────────────────────────────────────┐
│ User has:                          │
│  • cropseq_data.h5ad (from scanpy) │
│  • guides.csv (from lab)           │
└────────────────────────────────────┘
           │
           ▼
Step 3: Run analysis
┌────────────────────────────────────────────┐
│ $ perturbio analyze \                      │
│     cropseq_data.h5ad \                    │
│     --guides guides.csv                    │
│                                            │
│ [Progress output showing each step...]     │
│ ✓ Analysis Complete!                       │
│ Results saved to: perturbio_results_.../   │
└────────────────────────────────────────────┘
           │
           ▼
Step 4: Review results
┌────────────────────────────────────────────┐
│ $ open results/figures/volcano_MYC.png     │
│                                            │
│ [Beautiful volcano plot appears]           │
│                                            │
│ 😮 "Wow, MYC is knocked down and all      │
│     its targets are affected!"             │
└────────────────────────────────────────────┘
           │
           ▼
Step 5: Share with advisor
┌────────────────────────────────────────────┐
│ "The CRISPR knockdown worked! Look at     │
│  these results from Perturbio."            │
│                                            │
│ ✓ Mission accomplished                     │
└────────────────────────────────────────────┘
```

### Journey 2: Experienced User (Custom Analysis)

```
┌─────────────────────────────────────────────────────────────────┐
│ Persona: Computational biologist, wants fine control            │
│ Goal: Integrate Perturbio into existing scanpy workflow         │
│ Time: 15 minutes                                                 │
└─────────────────────────────────────────────────────────────────┘

Step 1: Load and preprocess with scanpy
┌────────────────────────────────────────────┐
│ import scanpy as sc                        │
│ import perturbio as pt                     │
│                                            │
│ adata = sc.read_h5ad("data.h5ad")         │
│ sc.pp.filter_cells(adata, min_genes=200)  │
│ sc.pp.normalize_total(adata)              │
│ sc.pp.log1p(adata)                        │
└────────────────────────────────────────────┘
           │
           ▼
Step 2: Extract guides with Perturbio
┌────────────────────────────────────────────┐
│ pt.guides.extract(                         │
│     adata,                                 │
│     guide_file="guides.csv",               │
│     min_umis=5,  # Custom threshold        │
│     key_added="perturbation"               │
│ )                                          │
│                                            │
│ print(adata.obs['perturbation'].value_    │
│       counts())                            │
└────────────────────────────────────────────┘
           │
           ▼
Step 3: Custom dimensionality reduction
┌────────────────────────────────────────────┐
│ sc.pp.highly_variable_genes(adata)        │
│ sc.tl.pca(adata, n_comps=50)              │
│ sc.pp.neighbors(adata, n_neighbors=15)    │
│ sc.tl.umap(adata, min_dist=0.3)           │
│                                            │
│ # Visualize perturbations                  │
│ sc.pl.umap(adata, color='perturbation')   │
└────────────────────────────────────────────┘
           │
           ▼
Step 4: Differential expression
┌────────────────────────────────────────────┐
│ pt.tl.differential_expression(             │
│     adata,                                 │
│     groupby='perturbation',                │
│     control='non-targeting',               │
│     min_cells=20,  # Stricter threshold    │
│     key_added='my_de_results'              │
│ )                                          │
└────────────────────────────────────────────┘
           │
           ▼
Step 5: Custom visualization
┌────────────────────────────────────────────┐
│ # Volcano plot for specific perturbation   │
│ pt.pl.volcano(                             │
│     adata,                                 │
│     perturbation='BRCA1_guide1',           │
│     de_key='my_de_results',                │
│     fdr_threshold=0.01,                    │
│     label_top=15                           │
│ )                                          │
│                                            │
│ # Export for downstream analysis           │
│ de_results = adata.uns['my_de_results']    │
│ de_results.to_csv("custom_de.csv")         │
└────────────────────────────────────────────┘
```

---

## File Structure & Outputs

### Input File Requirements

```
Guide Library CSV Format:
┌──────────────────────────────────────────────────────────────┐
│ guide_id          | target_gene | guide_sequence             │
├──────────────────────────────────────────────────────────────┤
│ BRCA1_guide1      | BRCA1       | GCACTCAGGAAACAGCTATG       │
│ BRCA1_guide2      | BRCA1       | CTGAAGACTGCTCAGTGTAG       │
│ MYC_guide1        | MYC         | GTACTTGGTGAGGCCAGCGC       │
│ non-targeting_1   | control     | GTAGCGAACGTGTCCGGCGT       │
└──────────────────────────────────────────────────────────────┘

Minimum required columns:
  • guide_id: Unique identifier
  • target_gene: Gene being targeted (or "control")

Optional columns:
  • guide_sequence: RNA sequence
  • pam_sequence: PAM site
  • chromosome, start, end: Genomic coordinates
```

### Output Directory Structure

```
perturbio_results_20250126_143022/
│
├── summary.txt                          # Human-readable summary
│   └── Analysis overview, cell counts, top findings
│
├── data/
│   └── annotated_data.h5ad              # Input + guide annotations
│       └── Contains:
│           • Original count matrix
│           • adata.obs['perturbation']
│           • adata.obs['guide_identity']
│           • adata.obs['guide_umi_count']
│           • adata.uns['perturbio_de']
│
├── results/
│   ├── differential_expression.csv      # All genes, all perturbations
│   │   └── Columns: perturbation, gene, log_fc, pval, pval_adj, rank
│   │
│   ├── top_hits_summary.csv             # Top 50 genes per perturbation
│   │
│   ├── perturbations.csv                # Cell → guide assignments
│   │   └── Columns: cell_barcode, guide_identity, perturbation,
│   │                umi_count, confidence
│   │
│   └── guide_statistics.csv             # Per-guide summary stats
│       └── Columns: guide_id, target_gene, n_cells, n_sig_genes,
│                    top_gene, top_log_fc
│
└── figures/
    ├── perturbation_counts.png          # Bar chart of cells per guide
    ├── umap_perturbations.png           # UMAP colored by guide
    │
    ├── volcano/                         # Volcano plots
    │   ├── volcano_BRCA1_guide1.png
    │   ├── volcano_MYC_guide1.png
    │   ├── volcano_TP53_guide1.png
    │   └── ...
    │
    └── combined_volcano.pdf             # Multi-panel figure (optional)
```

### Example Output Files

#### summary.txt
```
╔═══════════════════════════════════════════════════════════════╗
║             Perturbio Analysis Summary                        ║
╚═══════════════════════════════════════════════════════════════╝

Analysis Date: 2025-01-26 14:30:22
Input File: cropseq_data.h5ad
Guide Library: guides.csv

──────────────────────────────────────────────────────────────────
DATASET OVERVIEW
──────────────────────────────────────────────────────────────────
Total Cells: 8,467
Total Genes: 23,451

──────────────────────────────────────────────────────────────────
GUIDE ASSIGNMENT
──────────────────────────────────────────────────────────────────
Cells with guides: 8,102 (95.7%)
Non-targeting controls: 289 (3.4%)
Unassigned cells: 76 (0.9%)

Guides tested: 24
Cells per guide (median): 337

──────────────────────────────────────────────────────────────────
DIFFERENTIAL EXPRESSION
──────────────────────────────────────────────────────────────────
Total significant genes: 3,847 (FDR < 0.05)
Average genes per perturbation: 160

Top 5 Perturbations by Effect Size:
1. MYC_guide1       - 567 sig genes  - MYC ↓4.2 (FDR=1.2e-78)
2. TP53_guide1      - 689 sig genes  - TP53 ↓3.1 (FDR=3.4e-56)
3. BRCA1_guide1     - 421 sig genes  - BRCA1 ↓3.8 (FDR=7.8e-45)
4. EGFR_guide2      - 234 sig genes  - EGFR ↓2.9 (FDR=2.1e-34)
5. KRAS_guide1      - 189 sig genes  - KRAS ↓2.3 (FDR=5.6e-29)

──────────────────────────────────────────────────────────────────
OUTPUT FILES
──────────────────────────────────────────────────────────────────
✓ Differential expression results
✓ Annotated dataset with perturbations
✓ 27 visualization plots
✓ Summary statistics

Total analysis time: 2m 34s
```

#### differential_expression.csv (excerpt)
```
perturbation,gene,log_fc,pval,pval_adj,rank
BRCA1_guide1,BRCA1,-3.82,1.23e-45,7.84e-43,1
BRCA1_guide1,RAD51,-2.14,3.45e-23,1.12e-21,2
BRCA1_guide1,BRIP1,-1.89,7.82e-18,2.34e-16,3
BRCA1_guide1,PALB2,-1.67,2.31e-15,5.67e-14,4
MYC_guide1,MYC,-4.21,1.15e-78,8.92e-76,1
MYC_guide1,CDK4,-2.98,4.56e-45,2.34e-43,2
...
```

---

## Visual Design Elements

### Progress Indicators

```
Loading data:
  ↻ Loading... (spinner while working)
  ✓ Complete (checkmark when done)
  ✗ Failed (X on error)
  ⚠ Warning (warning symbol)

Progress bar:
  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% (24/24)

  Or with description:
  Analyzing perturbations ━━━━━━━━━━━━━━━━━━━ 87% (21/24)
```

### Tables

```
Simple table:
┌──────────────┬─────────┐
│ Header 1     │ Header 2│
├──────────────┼─────────┤
│ Value 1      │ Value 2 │
└──────────────┴─────────┘

Bold table (for summaries):
╔══════════════╦═════════╗
║ Header 1     ║ Header 2║
╠══════════════╬═════════╣
║ Value 1      ║ Value 2 ║
╚══════════════╩═════════╝
```

### Color Scheme (Terminal)

```
When terminal colors are supported:

🟢 Success messages (green)
🔵 Info messages (blue)
🟡 Warnings (yellow)
🔴 Errors (red)
⚪ Regular output (white/default)

Example:
  ✓ Analysis complete (green)
  ⚠ Low cell count for KRAS_guide2 (yellow)
  ✗ File not found (red)
```

---

## API Consistency Principles

### 1. Naming Conventions
```
Functions: snake_case
  extract_guides()
  differential_expression()

Classes: PascalCase
  CropSeqAnalyzer
  DEResults

Constants: UPPER_SNAKE_CASE
  DEFAULT_FDR_THRESHOLD = 0.05
  MIN_CELLS_PER_PERTURBATION = 10
```

### 2. Parameter Patterns
```
Common across all functions:

Input data:
  - adata: AnnData object (modified in-place for scanpy-style)
  - data: str | AnnData (for high-level API)

Output keys:
  - key_added: str (where to store results in adata.obs/uns)

Thresholds:
  - min_cells: int (minimum cells for analysis)
  - min_umis: int (minimum UMIs for guide assignment)
  - fdr_threshold: float (significance threshold)

Control:
  - control: str | list (control perturbation labels)
  - control_label: str (how to identify controls)
```

### 3. Return Value Conventions
```
Scanpy-style (in-place modification):
  pt.guides.extract(adata, ...)  # Returns None, modifies adata
  pt.tl.differential_expression(adata, ...)  # Returns None

High-level API (returns results):
  analyzer.run()  # Returns self for chaining
  analyzer.top_hits(...)  # Returns DataFrame

Plotting (returns figure):
  pt.pl.volcano(...)  # Returns matplotlib Figure
```

---

## Error Messages Design

### Error Message Template

```
┌─────────────────────────────────────────────────────────────────┐
│ ✗ [Error Type]                                                  │
├─────────────────────────────────────────────────────────────────┤
│ [Clear description of what went wrong]                          │
│                                                                 │
│ Context:                                                        │
│  • [Relevant info about user's input]                           │
│  • [What was expected vs what was found]                        │
│                                                                 │
│ Possible causes:                                                │
│  ✗ [Most likely cause]                                          │
│  ✗ [Alternative cause]                                          │
│                                                                 │
│ Solutions:                                                      │
│  → [Specific actionable fix]                                    │
│  → [Alternative fix]                                            │
│                                                                 │
│ For help: [documentation link]                                  │
└─────────────────────────────────────────────────────────────────┘
```

### Common Error Scenarios

```
File Not Found:
┌─────────────────────────────────────────────────────────────────┐
│ ✗ File Not Found                                                │
├─────────────────────────────────────────────────────────────────┤
│ Could not locate: /path/to/data.h5ad                            │
│                                                                 │
│ → Check file path is correct                                    │
│ → Ensure file extension is .h5ad                                │
│ → Try using absolute path instead of relative                   │
└─────────────────────────────────────────────────────────────────┘

Invalid Data Format:
┌─────────────────────────────────────────────────────────────────┐
│ ✗ Invalid AnnData Format                                        │
├─────────────────────────────────────────────────────────────────┤
│ Expected count matrix in adata.X, but found None                │
│                                                                 │
│ Your file contains:                                             │
│  • adata.layers: ['raw', 'normalized']                          │
│  • adata.X: None                                                │
│                                                                 │
│ → Specify which layer to use: --layer raw                       │
│ → Or copy layer to X: adata.X = adata.layers['raw']            │
└─────────────────────────────────────────────────────────────────┘

No Guides Detected:
[See earlier example in Terminal Output section]
```

---

## Summary: Design Philosophy

### Core Principles

1. **Simplicity First**
   - One command to run everything: `perturbio analyze data.h5ad`
   - Sensible defaults that work for 90% of use cases
   - Progressive disclosure: simple → advanced

2. **Clear Communication**
   - Show progress at every step
   - Explain what's happening (not just "processing...")
   - Helpful errors with specific solutions

3. **Scanpy-Compatible**
   - Works seamlessly with AnnData objects
   - Follows scanpy conventions (in-place modification)
   - Can be dropped into existing workflows

4. **Beautiful Output**
   - Publication-quality plots by default
   - Well-formatted terminal output
   - Organized results directory

5. **Fail Gracefully**
   - Detect problems early with clear errors
   - Suggest fixes, not just report failure
   - Validate inputs before long computations

---

**The Goal**: Make Crop-Seq analysis feel effortless. From raw data to biological insights in under 5 minutes, with beautiful visualizations and clear results. A tool that "just works" but gives power users full control when needed.

```
┌─────────────────────────────────────────────────────────────────┐
│                                                                 │
│  "The best interface is the one that gets out of your way."    │
│                                                                 │
│  Load data → Extract guides → Discover perturbations → Done.   │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```
