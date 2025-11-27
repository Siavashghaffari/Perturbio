# Perturbio Examples and Tutorials

This directory contains step-by-step tutorials and examples for using Perturbio to analyze Crop-Seq data.

## Tutorials

### 01. Quickstart Tutorial (5 minutes)
**File:** `01_quickstart.ipynb`
**Level:** Beginner
**Time:** 5-10 minutes

Learn the basics of Perturbio with a simple example:
- Load Crop-Seq data
- Extract guide barcodes
- Run differential expression
- Generate basic visualizations

Perfect for first-time users!

### 02. Complete Workflow Tutorial
**File:** `02_complete_workflow.ipynb`
**Level:** Intermediate
**Time:** 20-30 minutes

Comprehensive end-to-end analysis:
- Data preprocessing and QC
- Guide extraction with parameter tuning
- Differential expression analysis
- Advanced visualizations
- Results interpretation
- Exporting data and figures

### 03. Advanced Usage and Customization
**File:** `03_advanced_usage.ipynb`
**Level:** Advanced
**Time:** 30-45 minutes

Deep dive into advanced features:
- Scanpy integration
- Custom analysis workflows
- Parameter optimization
- Handling edge cases
- Batch processing multiple datasets
- Publication-quality figure customization

### 04. CLI Usage Examples
**File:** `04_cli_examples.md`
**Level:** All levels
**Time:** 10 minutes

Command-line interface examples:
- Basic CLI commands
- Batch processing scripts
- Integration with shell pipelines
- Automation workflows

## Example Data

The `data/` folder contains:
- `example_guide_library.csv` - Template guide library
- Instructions for obtaining example Crop-Seq datasets
- Small synthetic dataset for testing

## Running the Tutorials

### Option 1: Jupyter Notebook
```bash
cd examples
jupyter notebook
```

### Option 2: JupyterLab
```bash
cd examples
jupyter lab
```

### Option 3: VS Code
Open the `.ipynb` files directly in VS Code with the Jupyter extension.

## Prerequisites

Make sure you have Perturbio installed:
```bash
pip install perturbio
```

For the tutorials, you'll also need:
```bash
pip install jupyter matplotlib seaborn
```

## Getting Help

- 📖 [Full Documentation](https://perturbio.readthedocs.io)
- 💬 [GitHub Discussions](https://github.com/perturbio/perturbio/discussions)
- 🐛 [Report Issues](https://github.com/perturbio/perturbio/issues)

## Tutorial Structure

Each tutorial follows this format:
1. **Learning Objectives** - What you'll learn
2. **Setup** - Required imports and data
3. **Step-by-Step Analysis** - Detailed walkthrough with explanations
4. **Results Interpretation** - Understanding the outputs
5. **Next Steps** - What to explore next

Happy analyzing! 🧬
