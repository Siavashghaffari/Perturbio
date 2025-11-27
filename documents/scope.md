# Perturbio - Scope Document

## Overview

**Perturbio** is a comprehensive Python package for end-to-end analysis of Crop-Seq experiments (CRISPR pooled screens combined with single-cell RNA sequencing). The package enables researchers to process raw sequencing data, identify guide RNA perturbations in individual cells, quantify transcriptional responses, and discover genes that regulate cellular processes through systematic CRISPR perturbations at single-cell resolution.

**Target Users**: Computational biologists, bioinformaticians, and researchers conducting pooled CRISPR screens with single-cell readouts.

**Installation**: `pip install perturbio`

---

## Problem Statement

Crop-Seq experiments generate complex multi-modal data requiring:
- Identification of which guide RNA perturbed each cell (barcode extraction)
- Linking perturbations to transcriptional phenotypes (expression profiling)
- Statistical assessment of perturbation effects while controlling for technical/biological confounders
- Integration of multiple analysis modalities (differential expression, dimensionality reduction, pathway analysis)

Current workflows require stitching together disparate tools, custom scripts, and manual data wrangling. **Perturbio** provides a unified, standardized pipeline from raw data to biological insights.

---

## Core Features

### 1. Dataset Management
- **Public Data Access**: Download Crop-Seq datasets from GEO, SRA, and specialized repositories
- **Dataset Registry**: Curated catalog of published Crop-Seq experiments with metadata
- **Data Organization**: Automatic directory structure creation and file management
- **Version Control**: Track dataset versions and processing history

### 2. Preprocessing Pipeline
- **Quality Control Metrics**:
  - Per-cell QC (total counts, detected genes, mitochondrial percentage)
  - Per-guide QC (coverage, distribution across cells)
  - Doublet detection for cells with multiple guides
- **Read Alignment**: Interface to STAR, BWA, or Cell Ranger for genome alignment
- **Filtering**: Automated and customizable thresholds for low-quality cells/genes
- **Normalization**: Library size normalization, log-transformation, variance stabilization

### 3. CRISPR Barcode Extraction
- **Guide Detection**: Extract guide RNA sequences from sequencing reads (from transcriptome or dedicated guide library)
- **UMI Handling**: Count unique molecular identifiers for accurate guide quantification
- **Multi-guide Assignment**: Handle cells with multiple detected guides (MOI analysis)
- **Control Cell Identification**: Classify non-targeting control cells
- **Confidence Scoring**: Quality metrics for guide-cell assignments

### 4. Gene Expression Analysis
- **Differential Expression**:
  - Compare perturbed cells vs. controls (or other perturbations)
  - Handle batch effects and confounders
  - Support for multiple DE methods (Wilcoxon, t-test, DESeq2-like models)
- **Dose-Response Analysis**: Correlate guide abundance with expression changes
- **Target Gene Validation**: Verify knockdown/knockout efficiency for targeted genes
- **Expression Signatures**: Generate perturbation-specific gene signatures

### 5. Cell Cycle Scoring
- **Phase Assignment**: Classify cells into G1, S, G2/M phases
- **Species Support**: Built-in marker genes for human, mouse, rat
- **Custom Markers**: Allow user-defined cell cycle gene sets
- **Regression**: Option to regress out cell cycle effects from downstream analysis
- **Perturbation-Cycle Interaction**: Identify guides affecting cell proliferation

### 6. Dimensionality Reduction & Clustering
- **Embedding Generation**: UMAP, t-SNE, PCA, diffusion maps
- **Perturbation-Aware Embeddings**: Visualize cells colored by guide identity
- **Clustering**: Louvain, Leiden, hierarchical clustering
- **Trajectory Inference**: Pseudotime analysis for perturbations affecting differentiation
- **Batch Correction**: Harmony, Scanorama, or Combat integration

### 7. Statistical Testing
- **Perturbation Significance**:
  - Fisher exact test for categorical phenotypes
  - Permutation tests for expression differences
  - Energy distance for distributional changes
  - Mixscape-style classifier for perturbation detection
- **Multiple Testing Correction**: FDR (Benjamini-Hochberg), Bonferroni
- **Effect Size Estimation**: Cohen's d, log fold-change, percent change
- **Power Analysis**: Estimate sample sizes needed for future experiments

### 8. Gene Set Enrichment Analysis
- **Pathway Databases**: Integration with MSigDB, GO, KEGG, Reactome
- **Enrichment Methods**:
  - Over-representation analysis (ORA)
  - Gene Set Enrichment Analysis (GSEA)
  - Single-cell pathway scoring
- **Custom Gene Sets**: Support user-defined signatures
- **Visualization**: Enrichment plots, dot plots, network diagrams

### 9. Phenotype Classification
- **Supervised Learning**:
  - Train classifiers to predict cellular states from expression
  - Support for Random Forest, Logistic Regression, SVM, Neural Networks
- **Feature Selection**: Identify discriminative genes for phenotypes
- **Cross-Validation**: Robust performance evaluation
- **Perturbation Prediction**: Predict guide identity from expression profiles
- **Transfer Learning**: Apply trained models to new datasets

### 10. Visualization & Reporting
- **Automated Plots**:
  - QC violin plots and scatter plots
  - UMAP/t-SNE embeddings with perturbation overlays
  - Volcano plots for differential expression
  - Heatmaps of top perturbed genes
  - Guide coverage and distribution plots
  - Enrichment dot plots and bar charts
- **Publication Quality**: Customizable themes, export to PDF/PNG/SVG
- **Interactive Notebooks**: Auto-generated Jupyter notebooks with analysis steps
- **HTML Reports**: Standalone reports with interactive Plotly visualizations
- **Export Utilities**: Save figures in standardized directory structures

---

## Package Architecture

### Module Organization

```
perturbio/
├── data/              # Dataset download and management
├── preprocessing/     # QC, filtering, normalization
├── guides/            # CRISPR barcode extraction and assignment
├── expression/        # Differential expression analysis
├── cell_cycle/        # Cell cycle scoring and regression
├── reduction/         # Dimensionality reduction and clustering
├── statistics/        # Statistical testing for perturbations
├── enrichment/        # Gene set enrichment analysis
├── classification/    # ML models for phenotype prediction
├── visualization/     # Plotting utilities
├── reports/           # Automated report generation
├── io/                # Import/export utilities
├── utils/             # Helper functions
└── cli/               # Command-line interface
```

### Design Principles

1. **Modular**: Each component can be used independently or as part of pipeline
2. **Scanpy-Native**: All operations work with AnnData objects
3. **Reproducible**: Automatic logging of parameters and random seeds
4. **Extensible**: Plugin architecture for custom methods
5. **Performance**: Optimized for large datasets with sparse matrices and parallel processing

---

## User Interfaces

### 1. Python API
```python
import perturbio as pt

# Load data
adata = pt.read_10x_h5("data.h5")

# Preprocess
pt.pp.quality_control(adata)
pt.pp.normalize(adata)

# Extract guides
pt.guides.extract_from_transcriptome(adata)

# Analysis
pt.tl.differential_expression(adata, groupby='guide_identity')
pt.tl.cell_cycle_score(adata, species='human')
pt.tl.perturbation_significance(adata)

# Visualization
pt.pl.umap(adata, color='guide_identity')
pt.pl.volcano(adata, perturbation='BRCA1_guide1')

# Report
pt.reports.generate_html(adata, output="analysis_report.html")
```

### 2. Command-Line Interface
```bash
# Full pipeline
perturbio run \
  --input data.h5 \
  --guides guides.csv \
  --species human \
  --output results/

# Individual steps
perturbio preprocess --input data.h5 --output qc_data.h5
perturbio extract-guides --input qc_data.h5 --guides guides.csv
perturbio diff-exp --input processed.h5 --group guide_identity
perturbio report --input analyzed.h5 --output report.html
```

### 3. Configuration Files
- YAML/TOML configs for reproducible pipeline runs
- Parameter presets for common experimental designs

---

## Technical Requirements

### Core Dependencies
- **Python**: 3.9+
- **Data Processing**:
  - `scanpy` (>=1.9): Single-cell analysis foundation
  - `anndata` (>=0.8): Data structure
  - `pandas` (>=1.5): Tabular data
  - `numpy` (>=1.23): Numerical operations
  - `scipy` (>=1.9): Statistical functions
- **Machine Learning**:
  - `scikit-learn` (>=1.2): ML models and metrics
- **Visualization**:
  - `matplotlib` (>=3.6): Core plotting
  - `seaborn` (>=0.12): Statistical visualization
  - `plotly` (>=5.0): Interactive plots
- **Bioinformatics**:
  - `pysam` (>=0.20): BAM/SAM file processing
  - `pybedtools`: Genomic interval operations
- **CLI**:
  - `click` (>=8.0): Command-line interface
- **Utilities**:
  - `tqdm`: Progress bars
  - `joblib`: Parallel processing
  - `requests`: Dataset downloads

### Optional Dependencies
- `harmony-pytorch`: Batch correction
- `bbknn`: Graph-based batch correction
- `torch`: Deep learning models (if needed)
- `jupyter`: Notebook generation

### Performance Targets
- Handle 100,000+ cells without memory issues
- Utilize sparse matrices for gene expression
- Parallel processing for embarrassingly parallel tasks (per-guide DE)
- Leverage GPU for compatible operations (optional)

### Species Support
- **Built-in**: Human (hg38), Mouse (mm10), Rat (rn6)
- **Custom**: User-provided genome annotations (GTF/GFF)

### Experimental Designs
- **Pooled Screens**: Multiple guides in same culture
- **Arrayed Screens**: One guide per well (with multiplexing)
- **Single vs Multi-guide**: MOI analysis for cells with multiple perturbations
- **Time-Series**: Longitudinal Crop-Seq experiments

---

## Data Formats and Integration

### Input Formats
- **Raw Reads**: FASTQ (for guide extraction pipeline)
- **Count Matrices**:
  - 10x Genomics HDF5 (.h5)
  - Matrix Market (.mtx)
  - Loom (.loom)
  - Text (CSV/TSV)
- **Guide Annotations**: CSV/TSV with guide sequences and target genes
- **Cell Metadata**: CSV/TSV with experimental conditions, batches

### Output Formats
- **Primary**: AnnData HDF5 (.h5ad) - preserves all analysis layers
- **Expression Matrices**: CSV, TSV, Excel
- **Results Tables**:
  - Differential expression: CSV with statistics
  - Enrichment results: CSV with pathways and p-values
  - Cell metadata: CSV with assignments and scores
- **Figures**: PDF, PNG, SVG
- **Reports**: HTML with embedded interactive plots
- **Serialized Models**: Pickle/joblib for ML classifiers

### Integration Points
- **Scanpy Ecosystem**: Seamless interoperability with scanpy workflows
- **R Interoperability**: Export to Seurat-compatible formats (RDS via anndata2ri)
- **Track Hubs**: Generate UCSC/IGV tracks for guide locations
- **External Tools**: Interface with CellRanger, STAR, BWA outputs

---

## Deliverables

### 1. Core Package
- Well-documented codebase with type hints
- Comprehensive unit tests (pytest)
- Continuous integration (GitHub Actions)
- Code coverage >80%

### 2. Documentation
- **Installation Guide**: pip, conda, Docker options
- **API Reference**: Auto-generated from docstrings (Sphinx)
- **User Guide**: Conceptual overview of Crop-Seq analysis
- **Tutorials**:
  - Quickstart (30 min)
  - End-to-end analysis (2-3 hours)
  - Advanced workflows (custom analyses)
- **Example Datasets**: 3-5 small Crop-Seq datasets (<5K cells) for testing
- **Best Practices**: Recommendations for experimental design and quality control

### 3. Example Notebooks
- **Basic Analysis**: Load data → preprocess → extract guides → visualize
- **Differential Expression**: Find perturbed genes, pathway enrichment
- **Machine Learning**: Train classifier to predict perturbations
- **Comparative Analysis**: Multi-guide effects, interaction screens
- **Custom Workflows**: Extend Perturbio with user functions

### 4. Community Resources
- GitHub repository with issue tracker
- Example gallery (rendered notebooks)
- FAQ and troubleshooting guide
- Contributing guidelines for developers

---

## Success Criteria

### Functionality
- ✅ Successfully analyze published Crop-Seq datasets end-to-end
- ✅ Reproduce key findings from landmark papers (e.g., Dixit et al. 2016, Datlinger et al. 2017)
- ✅ Generate publication-quality figures with <10 lines of code
- ✅ Complete analysis of 50K cell dataset in <30 minutes on standard laptop

### Usability
- ✅ Users with scanpy experience can start analyzing within 1 hour
- ✅ CLI enables non-programmers to run standard pipelines
- ✅ Clear error messages with actionable suggestions
- ✅ Comprehensive logging for debugging

### Adoption
- ✅ 5+ external research groups adopt for real analyses within 6 months
- ✅ Cited in peer-reviewed publications
- ✅ Active community engagement (GitHub stars, issues, pull requests)

---

## Out of Scope (Initial Release)

### V1 Exclusions
- ❌ Web-based GUI (command-line and API only)
- ❌ Cloud-hosted analysis platform
- ❌ Real-time experiment monitoring
- ❌ Automated guide library design (focus on analysis, not design)
- ❌ Base editing or prime editing screens (standard CRISPR only)
- ❌ Spatial transcriptomics integration (scRNA-seq only)
- ❌ Multi-omics (CITE-seq, ATAC-seq) - may add in future

### Future Considerations
- **V2**: Optical pooled screens (image-based phenotyping)
- **V3**: Combinatorial perturbations (pairwise guide effects)
- **V4**: Integration with genetic interaction databases (STRING, BioGRID)

---

## Timeline and Milestones

### Phase 1: Core Infrastructure (Months 1-2)
- Dataset management and I/O
- Preprocessing pipeline
- Guide extraction

### Phase 2: Analysis Methods (Months 3-4)
- Differential expression
- Statistical testing
- Dimensionality reduction

### Phase 3: Advanced Features (Months 5-6)
- Cell cycle scoring
- Enrichment analysis
- ML classification

### Phase 4: Visualization & Reporting (Month 7)
- Plotting functions
- Automated reports
- Interactive notebooks

### Phase 5: Documentation & Release (Month 8)
- API documentation
- Tutorials and examples
- PyPI release

---

## Related Tools and Positioning

### Existing Tools
- **Mixscape** (Seurat): Perturbation classification, limited to R ecosystem
- **SCEPTRE**: Statistical framework, narrow focus on QTL mapping
- **Perturb-seq analysis scripts**: Custom pipelines in individual papers

### Perturbio Advantages
- **Unified Pipeline**: End-to-end workflow in single package
- **Python/Scanpy Native**: Integrates with dominant single-cell ecosystem
- **Flexible**: Both CLI and programmatic access
- **Comprehensive**: Covers full analysis spectrum (QC → enrichment → reporting)
- **Reproducible**: Automated logging and configuration management
- **Scalable**: Optimized for large datasets

---

## Development Philosophy

1. **User-Centric**: Design API based on actual Crop-Seq analysis workflows
2. **Defaults That Work**: Sensible parameter defaults for 80% of use cases
3. **Transparency**: Users should understand what each function does
4. **Reproducibility First**: Track all parameters and random seeds automatically
5. **Community-Driven**: Accept contributions, respond to issues promptly
6. **Quality Over Features**: Robust core > sprawling half-baked features

---

## Summary

Perturbio aims to be the definitive Python package for Crop-Seq analysis, providing researchers with a powerful, user-friendly toolkit to extract biological insights from CRISPR-based single-cell screens. By integrating seamlessly with the scanpy ecosystem and offering both programmatic and command-line interfaces, Perturbio will accelerate discoveries in functional genomics and systems biology.

**Key Differentiator**: The only comprehensive, production-ready, Python-native solution for end-to-end Crop-Seq analysis.