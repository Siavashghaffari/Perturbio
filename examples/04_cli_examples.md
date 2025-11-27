# CLI Usage Examples

**Time:** 10 minutes
**Level:** All levels

## Overview

Perturbio provides a powerful command-line interface for batch processing, automation, and integration into computational pipelines. This tutorial shows you how to use Perturbio from the terminal.

---

## Basic Commands

### Get Help

```bash
# View all available commands
perturbio --help

# Get help for a specific command
perturbio analyze --help

# Check installed version
perturbio --version
```

### Simple Analysis

```bash
# Basic analysis with default parameters
perturbio analyze cropseq_data.h5ad \
  --guides guide_library.csv \
  --output results/

# Specify control label
perturbio analyze cropseq_data.h5ad \
  --guides guide_library.csv \
  --control non-targeting \
  --output experiment_001/
```

---

## Advanced Options

### Custom Parameters

```bash
# Fine-tune analysis parameters
perturbio analyze cropseq_data.h5ad \
  --guides guide_library.csv \
  --output results/ \
  --control non-targeting \
  --min-cells 20 \
  --min-umis 5 \
  --fdr 0.01
```

**Parameters:**
- `--control`: Control perturbation label (default: 'non-targeting')
- `--min-cells`: Minimum cells per perturbation (default: 10)
- `--min-umis`: Minimum UMIs to assign a guide (default: 3)
- `--fdr`: FDR threshold for significance (default: 0.05)

---

## Batch Processing

### Process Multiple Datasets

Create a bash script to analyze multiple experiments:

```bash
#!/bin/bash
# batch_analyze.sh

# List of datasets
DATASETS=(
  "experiment_001.h5ad"
  "experiment_002.h5ad"
  "experiment_003.h5ad"
)

GUIDE_LIB="shared_guide_library.csv"

# Process each dataset
for dataset in "${DATASETS[@]}"; do
  echo "Processing $dataset..."

  # Extract experiment name
  exp_name=$(basename "$dataset" .h5ad)

  # Run analysis
  perturbio analyze "$dataset" \
    --guides "$GUIDE_LIB" \
    --output "results_${exp_name}/" \
    --control non-targeting \
    --min-cells 15 \
    --fdr 0.05

  echo "✓ Completed $dataset"
done

echo "All datasets processed!"
```

Run the script:
```bash
chmod +x batch_analyze.sh
./batch_analyze.sh
```

---

## Integration with Shell Pipelines

### Pre-process → Analyze → Summarize

```bash
#!/bin/bash
# full_pipeline.sh

INPUT="raw_cropseq.h5ad"
GUIDES="guides.csv"
OUTPUT="final_results"

# Step 1: Quality control with scanpy (Python)
python - <<EOF
import scanpy as sc
adata = sc.read_h5ad("$INPUT")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.write("qc_data.h5ad")
print("✓ QC complete")
EOF

# Step 2: Perturbio analysis
echo "Running Perturbio analysis..."
perturbio analyze qc_data.h5ad \
  --guides "$GUIDES" \
  --output "$OUTPUT" \
  --min-cells 20 \
  --fdr 0.05

# Step 3: Summary statistics
echo "Generating summary..."
python - <<EOF
import pandas as pd
de = pd.read_csv("$OUTPUT/differential_expression.csv")
summary = de.groupby('perturbation').agg({
    'gene': 'count',
    'pval_adj': lambda x: (x < 0.05).sum()
})
summary.columns = ['Total_Genes', 'Significant_Genes']
print("\nSummary Statistics:")
print(summary)
summary.to_csv("$OUTPUT/summary_stats.csv")
EOF

# Cleanup
rm qc_data.h5ad

echo "✓ Pipeline complete! Results in $OUTPUT/"
```

---

## Parallel Processing

### Analyze Multiple Experiments in Parallel

Using GNU Parallel:

```bash
#!/bin/bash
# parallel_analyze.sh

# Function to analyze one dataset
analyze_dataset() {
  local dataset=$1
  local guide_lib=$2
  local exp_name=$(basename "$dataset" .h5ad)

  echo "Starting $exp_name..."

  perturbio analyze "$dataset" \
    --guides "$guide_lib" \
    --output "results_${exp_name}/" \
    --min-cells 15

  echo "✓ Finished $exp_name"
}

# Export function for parallel
export -f analyze_dataset

# Run in parallel (4 jobs at a time)
ls experiments/*.h5ad | \
  parallel -j 4 analyze_dataset {} guide_library.csv

echo "All analyses complete!"
```

Install GNU Parallel:
```bash
# macOS
brew install parallel

# Linux (Ubuntu/Debian)
sudo apt-get install parallel
```

---

## Automation Workflows

### Automated Nightly Analysis

Create a cron job to run analysis automatically:

```bash
# Edit crontab
crontab -e

# Add this line to run daily at 2 AM
0 2 * * * /path/to/your/analyze_script.sh >> /path/to/logs/perturbio.log 2>&1
```

Example script for automated runs:

```bash
#!/bin/bash
# automated_analysis.sh

DATE=$(date +%Y%m%d)
LOG_FILE="logs/perturbio_${DATE}.log"

echo "Starting automated analysis: $(date)" | tee -a "$LOG_FILE"

# Check for new data
NEW_DATA=$(find data/incoming/ -name "*.h5ad" -mtime -1)

if [ -z "$NEW_DATA" ]; then
  echo "No new data found" | tee -a "$LOG_FILE"
  exit 0
fi

# Process new data
for dataset in $NEW_DATA; do
  echo "Processing $dataset" | tee -a "$LOG_FILE"

  exp_name=$(basename "$dataset" .h5ad)

  perturbio analyze "$dataset" \
    --guides guide_library.csv \
    --output "results_${DATE}_${exp_name}/" \
    2>&1 | tee -a "$LOG_FILE"

  # Move processed data
  mv "$dataset" "data/processed/"
done

echo "Analysis complete: $(date)" | tee -a "$LOG_FILE"

# Send notification (optional)
# echo "Perturbio analysis complete" | mail -s "Analysis Report" you@email.com
```

---

## Output Organization

### Organizing Results by Experiment

```bash
#!/bin/bash
# organize_results.sh

# Create organized directory structure
EXPERIMENT="CRISPRi_Screen_2024"
DATE=$(date +%Y%m%d)
RESULT_DIR="results/${EXPERIMENT}/${DATE}"

mkdir -p "$RESULT_DIR"/{raw,analysis,figures,reports}

# Run analysis
perturbio analyze data/experiment.h5ad \
  --guides guides.csv \
  --output "$RESULT_DIR/analysis/"

# Copy raw data
cp data/experiment.h5ad "$RESULT_DIR/raw/"
cp guides.csv "$RESULT_DIR/raw/"

# Move figures
mv "$RESULT_DIR/analysis/"*.png "$RESULT_DIR/figures/"

# Generate report
cat > "$RESULT_DIR/reports/README.md" <<EOF
# $EXPERIMENT Analysis Report
**Date:** $(date)
**Analyst:** $(whoami)

## Files
- Raw data: \`raw/experiment.h5ad\`
- Guide library: \`raw/guides.csv\`
- DE results: \`analysis/differential_expression.csv\`
- Figures: \`figures/*.png\`

## Parameters
- Min cells: 10
- Min UMIs: 3
- FDR threshold: 0.05
EOF

echo "✓ Results organized in $RESULT_DIR"
```

---

## Error Handling

### Robust Pipeline with Error Checking

```bash
#!/bin/bash
# robust_pipeline.sh

set -euo pipefail  # Exit on error, undefined vars, pipe failures

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

error_exit() {
  echo -e "${RED}ERROR: $1${NC}" >&2
  exit 1
}

success_msg() {
  echo -e "${GREEN}✓ $1${NC}"
}

# Check inputs exist
[ -f "data.h5ad" ] || error_exit "data.h5ad not found"
[ -f "guides.csv" ] || error_exit "guides.csv not found"

# Run analysis with error checking
echo "Starting analysis..."
if perturbio analyze data.h5ad \
     --guides guides.csv \
     --output results/ \
     --min-cells 10; then
  success_msg "Analysis complete"
else
  error_exit "Perturbio analysis failed"
fi

# Verify outputs
[ -f "results/differential_expression.csv" ] || \
  error_exit "Expected output file not found"

success_msg "All outputs generated successfully"
```

---

## Performance Monitoring

### Track Execution Time and Resources

```bash
#!/bin/bash
# monitor_analysis.sh

DATASET="large_cropseq.h5ad"
LOG="performance.log"

echo "=== Analysis Performance Report ===" | tee "$LOG"
echo "Dataset: $DATASET" | tee -a "$LOG"
echo "Start time: $(date)" | tee -a "$LOG"

# Time the analysis
/usr/bin/time -v perturbio analyze "$DATASET" \
  --guides guides.csv \
  --output results/ \
  2>&1 | tee -a "$LOG"

echo "End time: $(date)" | tee -a "$LOG"

# Extract key metrics (Linux)
if command -v grep &> /dev/null; then
  echo -e "\n=== Resource Usage ===" | tee -a "$LOG"
  grep "Maximum resident set size" "$LOG" || true
  grep "Elapsed (wall clock) time" "$LOG" || true
fi
```

---

## Integration with Workflow Managers

### Snakemake Workflow

Create `Snakefile`:

```python
# Snakefile for Perturbio pipeline

SAMPLES = ["exp001", "exp002", "exp003"]

rule all:
    input:
        expand("results/{sample}/differential_expression.csv", sample=SAMPLES)

rule analyze:
    input:
        data="data/{sample}.h5ad",
        guides="guides.csv"
    output:
        de="results/{sample}/differential_expression.csv",
        summary="results/{sample}/summary.txt"
    log:
        "logs/{sample}.log"
    shell:
        """
        perturbio analyze {input.data} \
          --guides {input.guides} \
          --output results/{wildcards.sample}/ \
          2> {log}
        """

rule summarize:
    input:
        expand("results/{sample}/differential_expression.csv", sample=SAMPLES)
    output:
        "final_summary.csv"
    run:
        import pandas as pd
        dfs = []
        for f in input:
            df = pd.read_csv(f)
            sample = f.split('/')[1]
            df['sample'] = sample
            dfs.append(df)
        pd.concat(dfs).to_csv(output[0], index=False)
```

Run with:
```bash
snakemake --cores 4
```

---

## Common Use Cases

### 1. Quick Screen Analysis

```bash
# Fast analysis with lenient parameters
perturbio analyze screen.h5ad \
  --guides guides.csv \
  --output quick_results/ \
  --min-cells 5 \
  --fdr 0.1
```

### 2. High-Confidence Hits

```bash
# Stringent parameters for publication
perturbio analyze screen.h5ad \
  --guides guides.csv \
  --output stringent_results/ \
  --min-cells 50 \
  --min-umis 5 \
  --fdr 0.01
```

### 3. Compare Parameters

```bash
#!/bin/bash
# compare_parameters.sh

# Test different min-cells thresholds
for min_cells in 10 20 30 50; do
  perturbio analyze data.h5ad \
    --guides guides.csv \
    --output "results_mincells_${min_cells}/" \
    --min-cells "$min_cells"
done

echo "Parameter sweep complete!"
```

---

## Debugging

### Verbose Output

```bash
# Python's verbose mode
python -v -m perturbio.cli analyze data.h5ad \
  --guides guides.csv \
  --output results/

# Or use Python directly for full traceback
python <<EOF
import perturbio as pt
analyzer = pt.CropSeqAnalyzer('data.h5ad')
results = analyzer.run('guides.csv')
analyzer.export('results/')
EOF
```

### Check Package Installation

```bash
# Verify installation
python -c "import perturbio; print(perturbio.__version__)"

# List installed files
pip show -f perturbio

# Check dependencies
pip list | grep -E "scanpy|anndata|pandas"
```

---

## Tips and Best Practices

### 1. **Use Absolute Paths** in Scripts
```bash
# Good
perturbio analyze /path/to/data.h5ad --output /path/to/results/

# Avoid relative paths in cron jobs
```

### 2. **Log Everything**
```bash
# Redirect output to log file
perturbio analyze data.h5ad \
  --guides guides.csv \
  --output results/ \
  > analysis.log 2>&1
```

### 3. **Validate Inputs First**
```bash
# Check file exists and is valid
python -c "import scanpy as sc; sc.read_h5ad('data.h5ad')" || exit 1
```

### 4. **Use Environment Variables**
```bash
export PERTURBIO_GUIDES="guide_library.csv"
export PERTURBIO_OUTPUT="results/"
export PERTURBIO_MIN_CELLS=20

perturbio analyze data.h5ad \
  --guides "$PERTURBIO_GUIDES" \
  --output "$PERTURBIO_OUTPUT" \
  --min-cells "$PERTURBIO_MIN_CELLS"
```

---

## Summary

You've learned how to:

✅ Run Perturbio from the command line
✅ Batch process multiple datasets
✅ Integrate into shell pipelines
✅ Automate analyses with cron
✅ Monitor performance
✅ Handle errors robustly
✅ Use workflow managers (Snakemake)

## Next Steps

- Check out **Tutorial 01** for Python API basics
- See **Tutorial 03** for advanced scanpy integration
- Read the [full documentation](https://perturbio.readthedocs.io)

Happy automating! 🚀
