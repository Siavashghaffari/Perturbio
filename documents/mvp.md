# Perturbio MVP - The Magical Moment

## Vision

**The Magic**: A researcher drops in their Crop-Seq dataset and within 5 minutes sees exactly which genes are significantly affected by each CRISPR perturbation, with beautiful visualizations that tell the story.

**Core Workflow**:
```
Raw Crop-Seq Data → Extract Guide Barcodes → Identify Perturbed Genes → Visualize Results
```

**MVP Mantra**: Skip the complexity. Deliver the insight.

---

## MVP Scope: The Essentials

### What We Build

#### 1. Data Input (2 formats)
- **H5AD files**: Pre-processed AnnData objects from scanpy/Cell Ranger
- **BAM files** (stretch): Aligned reads for guide extraction
- Single dataset at a time (no batch processing)

#### 2. Guide Barcode Extraction
- Read CRISPR guide sequences from cell barcodes
- Assign each cell to its perturbation (guide identity)
- Handle non-targeting controls
- Simple assignment: one dominant guide per cell (ignore multi-guide for now)
- Output: Add `guide_identity` and `perturbation` columns to `adata.obs`

#### 3. Differential Expression
- Compare each perturbation vs. non-targeting controls
- **Method**: Wilcoxon rank-sum test (scanpy's built-in)
- Compute log fold-change and adjusted p-values (Benjamini-Hochberg FDR)
- Output: DataFrame with genes ranked by significance for each perturbation

#### 4. Visualizations (The 3 Essential Plots)
1. **Perturbation Assignment Plot**: Bar chart showing number of cells per guide
2. **Volcano Plot**: Top hits for a selected perturbation (log FC vs -log10 p-value)
3. **UMAP with Perturbations**: Dimensionality reduction colored by guide identity

#### 5. Results Export
- **CSV**: Differential expression results for all perturbations
- **H5AD**: Original data + added perturbation annotations
- **Figures**: PNG/PDF of the 3 key plots

---

## User Experience

### Command-Line Interface

```bash
# Simplest possible usage
perturbio analyze cropseq_data.h5ad

# With options
perturbio analyze \
  --input cropseq_data.h5ad \
  --guides guides.csv \
  --output results/ \
  --control-label "non-targeting"

# Output structure:
# results/
#   ├── perturbations.csv           # Cell → guide assignments
#   ├── differential_expression.csv # All perturbations, all genes
#   ├── top_hits_summary.csv        # Top 50 genes per perturbation
#   ├── perturbation_counts.png     # Bar chart
#   ├── volcano_BRCA1.png           # Volcano for each perturbation
#   └── umap_perturbations.png      # UMAP colored by guide
```

### Python API

```python
from perturbio import CropSeqAnalyzer

# One-liner magic
analyzer = CropSeqAnalyzer("cropseq_data.h5ad")
results = analyzer.run()  # Does everything: extract → DE → plot

# Step-by-step control
analyzer = CropSeqAnalyzer("cropseq_data.h5ad")
analyzer.extract_guides(guide_file="guides.csv")
analyzer.differential_expression(control_label="non-targeting")
analyzer.plot_volcano(perturbation="BRCA1_guide1")
analyzer.plot_umap(color_by="perturbation")
analyzer.export(output_dir="results/")

# Access results
print(analyzer.top_hits("BRCA1_guide1", n=20))  # Top 20 perturbed genes
print(analyzer.cells_per_perturbation)  # Cell counts
```

### Scanpy Integration (Advanced Users)

```python
import scanpy as sc
import perturbio as pt

# Standard scanpy workflow
adata = sc.read_h5ad("cropseq_data.h5ad")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

# Add Perturbio step
pt.extract_guides(adata, guide_file="guides.csv")  # Modifies adata in-place

# Continue with scanpy
sc.tl.pca(adata)
sc.tl.umap(adata)
sc.pl.umap(adata, color='perturbation')  # Color by extracted perturbations

# Differential expression the Perturbio way
pt.tl.differential_expression(adata, groupby='perturbation', control='non-targeting')

# Results stored in adata.uns['perturbio_de']
```

---

## MVP Architecture

### Package Structure (Minimal)

```
perturbio/
├── __init__.py
├── io/
│   ├── __init__.py
│   └── readers.py          # Load h5ad, BAM (future)
├── guides/
│   ├── __init__.py
│   └── extraction.py       # Extract and assign guides
├── analysis/
│   ├── __init__.py
│   └── differential.py     # Wilcoxon DE test
├── plotting/
│   ├── __init__.py
│   └── core.py             # 3 essential plots
├── core.py                 # CropSeqAnalyzer class
└── cli.py                  # Click-based CLI

tests/
├── test_extraction.py
├── test_differential.py
└── test_integration.py

examples/
├── example_data/
│   └── demo_cropseq.h5ad   # Small test dataset (1000 cells)
└── quickstart.ipynb        # 5-minute tutorial
```

### Key Classes

```python
# Core analyzer class
class CropSeqAnalyzer:
    def __init__(self, data: str | AnnData)
    def extract_guides(self, guide_file: str = None) -> Self
    def differential_expression(self, control_label: str = "non-targeting") -> Self
    def run(self) -> Self  # Run full pipeline
    def plot_volcano(self, perturbation: str) -> Figure
    def plot_umap(self, color_by: str = "perturbation") -> Figure
    def plot_perturbation_counts() -> Figure
    def top_hits(self, perturbation: str, n: int = 50) -> pd.DataFrame
    def export(self, output_dir: str) -> None

# Results container
@dataclass
class DEResults:
    genes: pd.DataFrame          # All genes, all perturbations
    perturbations: List[str]
    control_label: str

    def top_hits(self, perturbation: str, n: int = 50) -> pd.DataFrame
    def significant_genes(self, fdr_threshold: float = 0.05) -> pd.DataFrame
```

---

## Technical Implementation

### Core Dependencies (Minimal)
- `anndata` (>=0.8): Data structure
- `scanpy` (>=1.9): Single-cell analysis, built-in DE
- `pandas` (>=1.5): Results tables
- `numpy` (>=1.23): Numerical ops
- `matplotlib` (>=3.6): Plotting
- `seaborn` (>=0.12): Pretty plots
- `click` (>=8.0): CLI
- `pysam` (>=0.20): BAM reading (optional for MVP)

### Guide Extraction Strategy

**Input**:
- Guide library CSV with columns: `guide_sequence`, `target_gene`, `guide_id`
- AnnData object with cells × genes

**Algorithm**:
1. Look for guide RNA names in `adata.var_names` (gene names)
2. For each cell, identify guides with non-zero counts
3. Assign cell to guide with highest UMI count (simple winner-takes-all)
4. If no guides detected → label as "unassigned"
5. If multiple guides with similar counts → flag as "multiplet" (exclude for now)

**Output**:
- `adata.obs['guide_identity']`: Guide ID assigned to each cell
- `adata.obs['perturbation']`: Target gene for that guide
- `adata.obs['guide_umi_count']`: UMI count for assigned guide
- `adata.uns['guide_library']`: Guide metadata

### Differential Expression Strategy

**Method**: Scanpy's `sc.tl.rank_genes_groups()` with Wilcoxon
- For each perturbation: compare perturbed cells vs. non-targeting control cells
- Compute log fold-change, p-value, adjusted p-value (FDR)
- Rank genes by statistical significance

**Output**:
- DataFrame with columns: `perturbation`, `gene`, `log_fc`, `pval`, `pval_adj`, `rank`

---

## Success Criteria (MVP)

### Functional Requirements
✅ **Under 5 Minutes**: Full analysis (1000-10000 cells) completes in <5 min on laptop
✅ **Accurate Guide Assignment**: >95% agreement with manual inspection on test dataset
✅ **Reproducible DE Results**: Match scanpy's output when run manually
✅ **Publication-Quality Plots**: Figures ready for papers with minimal tweaking

### User Experience
✅ **Zero-Config Start**: `perturbio analyze data.h5ad` works out of the box
✅ **Clear Errors**: If guide assignment fails, tell user exactly why
✅ **Intuitive API**: Users familiar with scanpy can use Perturbio in <30 min

### Code Quality
✅ **Tests Pass**: 80%+ code coverage
✅ **Documented**: Every public function has docstring with example
✅ **Typed**: Type hints for all public APIs

---

## What's NOT in MVP

### Features to Skip (For Now)
❌ **Cell Cycle Scoring**: Assume confounders are minimal or pre-regressed
❌ **Advanced Statistics**: No Fisher exact, energy distance, permutation tests
❌ **Gene Set Enrichment**: No pathway analysis (just gene lists)
❌ **ML Classification**: No phenotype classifiers or perturbation predictors
❌ **Batch Correction**: Single dataset only, no integration
❌ **Multi-Guide Analysis**: Cells with multiple guides are flagged but not analyzed
❌ **Interactive Reports**: Just static plots (no HTML dashboards)
❌ **Dataset Management**: No download utilities for public data
❌ **Custom Genomes**: Human/mouse only, use pre-defined guide libraries

### Simplifications
- **One statistical test**: Wilcoxon only (scanpy's default)
- **One-vs-control design**: No perturbation-vs-perturbation comparisons
- **Simple guide assignment**: Highest UMI count wins
- **Pre-processed data**: Assume QC and normalization done by user

---

## MVP Development Timeline

### Week 1-2: Core Infrastructure
- [ ] Project setup (repo, dependencies, CI)
- [ ] AnnData I/O utilities
- [ ] Guide extraction algorithm
- [ ] Unit tests for guide assignment

### Week 3-4: Analysis Pipeline
- [ ] Differential expression wrapper
- [ ] Results data structures
- [ ] Integration tests with real Crop-Seq data

### Week 5-6: Visualization & CLI
- [ ] 3 core plotting functions
- [ ] CLI with Click
- [ ] CropSeqAnalyzer high-level API

### Week 7: Documentation & Polish
- [ ] API documentation
- [ ] Quickstart notebook
- [ ] Example dataset (1000 cells)
- [ ] README with installation + usage

### Week 8: Beta Testing
- [ ] External user testing (2-3 researchers)
- [ ] Bug fixes based on feedback
- [ ] Performance optimization

**Target MVP Release**: 8 weeks

---

## Example Workflow (The Magic in Action)

### Scenario: User has Crop-Seq data from Perturb-seq paper

```bash
# Download example data (pre-processed h5ad from published study)
wget https://example.com/dixit2016_cropseq.h5ad
wget https://example.com/dixit2016_guides.csv

# Run Perturbio
perturbio analyze dixit2016_cropseq.h5ad --guides dixit2016_guides.csv

# Output (after ~3 minutes):
# ✓ Loaded 7,738 cells × 32,738 genes
# ✓ Extracted guide identities for 7,420 cells (96%)
# ✓ Identified 24 perturbations + 1 control
# ✓ Differential expression: 24 perturbations vs non-targeting
# ✓ Found 1,847 significantly perturbed genes (FDR < 0.05)
# ✓ Generated 3 plots
# ✓ Results saved to results/
#
# Top perturbed gene: CXCR4_guide1 → CXCR4 (log2FC = -3.2, FDR = 1e-45) ✓
```

**User opens `results/volcano_CXCR4_guide1.png`**: Beautiful volcano plot showing CXCR4 knocked down as expected, plus downstream targets.

**The Magical Moment**: "It works! I can see my perturbations working!"

---

## Open Questions for MVP

### Technical Decisions
1. **Guide detection threshold**: How many UMIs needed to confidently assign a guide?
   - *Proposal*: Require ≥3 UMIs, flag cells with <3 as "low_confidence"

2. **Multiple guides per cell**: What's the MOI assumption?
   - *Proposal*: If cell has 2+ guides with >20% of max UMI count → label "multiplet", exclude from DE

3. **Minimum cells per perturbation**: Need enough for statistical power
   - *Proposal*: Require ≥10 cells, warn if <30 cells

4. **Control label**: How to identify non-targeting cells?
   - *Proposal*: Look for guides with "non-targeting", "control", "NTC" in name, or accept `--control-label` flag

### User Experience
1. **Default output location**: Current directory or subfolder?
   - *Proposal*: Create `perturbio_results_YYYYMMDD_HHMMSS/` to avoid overwriting

2. **Verbosity**: How much logging?
   - *Proposal*: Progress bars for long steps, summary stats at end, `--verbose` flag for details

3. **Error handling**: What if no guides detected?
   - *Proposal*: Check if any `adata.var_names` match guide library, show clear error with suggestions

---

## Post-MVP Roadmap (V0.2 - V1.0)

Once MVP proves the concept, add:

### V0.2 (Weeks 9-12)
- Cell cycle scoring and regression
- Batch processing (analyze multiple datasets)
- Perturbation-vs-perturbation comparisons

### V0.3 (Weeks 13-16)
- Gene set enrichment analysis (MSigDB integration)
- Additional statistical tests (Fisher exact, permutation)
- HTML reports with interactive plots

### V0.4 (Weeks 17-20)
- ML-based perturbation classification (Mixscape-style)
- Multi-guide analysis for MOI >1 experiments
- Performance optimization for 100K+ cells

### V1.0 (Weeks 21-24)
- Dataset management (download from GEO/SRA)
- Custom genome support
- Comprehensive documentation and tutorials
- PyPI release

---

## MVP Mantra (Repeat Daily)

**Focus**: One dataset, one analysis, three plots.

**Speed**: 5 minutes from data to insight.

**Magic**: The moment a researcher sees their perturbations working.

**Simplicity**: If it's not essential for the magic, it's not in the MVP.

---

## Success Looks Like

A graduate student with a Crop-Seq dataset:
1. Installs: `pip install perturbio` (30 seconds)
2. Runs: `perturbio analyze my_data.h5ad` (3 minutes)
3. Opens: `results/volcano_MYC_guide1.png`
4. Reacts: "Holy shit, it worked! MYC is knocked down and all the expected targets are affected!"
5. Shares: Posts figure on lab Slack with "Perturbio is amazing, check this out"

That's the MVP.