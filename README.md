# Proteomics MaxQuant LFQ Analysis Pipeline

This repository provides an end-to-end analysis workflow for MaxQuant label-free proteomics exported as `proteinGroups.txt`. It performs standard MaxQuant filtering, builds an LFQ intensity matrix using an explicit sample sheet, applies normalization and missing-value handling, runs differential abundance testing with limma, and generates core QC and results plots. Optional GO enrichment is included.

## Contents
- Repository structure
- Requirements
- Installation and reproducibility
- Inputs
  - proteinGroups.txt
  - metadata/sample_sheet.csv
  - metadata/contrasts.csv
- Run the pipeline
- Outputs
- Methods
- Troubleshooting
- License

---

## Repository structure

proteomics-maxquant-lfq/
  README.md
  run_pipeline.R
  renv.lock
  .gitignore
  R/
    utils.R
    01_import_qc.R
    02_normalize_impute.R
    03_differential_limma.R
    04_plots.R
    05_enrichment.R
  metadata/
    sample_sheet.csv
    contrasts.csv
  data_raw/
    proteinGroups.txt
  results/
    differential/
    enrichment/
  figures/

---

## Requirements
- R >= 4.2 recommended
- MaxQuant output includes LFQ intensity columns (commonly `LFQ intensity <sample>`)

---

## Installation and reproducibility

This project uses `renv` to manage package versions.

In R:

install.packages("renv")
renv::restore()

If you are creating this repository from scratch (first time only):

install.packages("renv")
renv::init()
# Install required packages used in scripts under R/
renv::snapshot()

---

## Inputs

### proteinGroups.txt
Place your MaxQuant file here:

data_raw/proteinGroups.txt

Notes:
- Avoid committing large raw files to GitHub. Keep `data_raw/` in `.gitignore`.

### metadata/sample_sheet.csv
This file explicitly maps each biological sample to the exact LFQ intensity column in `proteinGroups.txt` and defines the grouping used for modeling.

Required columns:
- sample: a short sample name used in outputs (must be unique)
- condition: group label used in limma design (for example Control, Treatment)
- batch: batch/run identifier (recommended even if all 1)
- mq_column: the exact column header from `proteinGroups.txt` for that sample’s LFQ values

Example:

sample,condition,batch,mq_column
N1,Normal,1,LFQ intensity N1
N2,Normal,1,LFQ intensity N2
N3,Normal,1,LFQ intensity N3
N4,Normal,1,LFQ intensity N4
S1,STIL,1,LFQ intensity S1
S2,STIL,1,LFQ intensity S2
S3,STIL,1,LFQ intensity S3
S4,STIL,1,LFQ intensity S4

Important:
- mq_column must match the column name in proteinGroups.txt exactly (including spaces/case).
- Column order does not matter because mapping is explicit.
- If your MaxQuant columns use a different prefix (for example Intensity instead of LFQ intensity), set mq_column accordingly.

### metadata/contrasts.csv
Each row defines a comparison to test. Conditions must match exactly the condition values used in sample_sheet.csv.

Columns:
- contrast: output label
- numerator: condition name for group 1
- denominator: condition name for group 2

Example:

contrast,numerator,denominator
Normal_vs_STIL,Normal,STIL
Normal_vs_STIC,Normal,STIC
Normal_vs_p53,Normal,p53

---

## Run the pipeline

From the repository root in R:

source("run_pipeline.R")

The pipeline will create (if missing):
- results/
- figures/

---

## Outputs

### Differential abundance results (limma)
Saved to:

results/differential/<contrast>.csv

Each file is a limma results table (logFC, P.Value, adj.P.Val, etc.) for the named contrast.

### Figures
Saved to figures/:
- PCA.png
- sample_correlation.png
- volcano_<contrast>.png
- heatmap_top50_<contrast>.png

### Optional enrichment
Saved to results/enrichment/:
- GO_BP_<contrast>.png
- GO_BP_<contrast>.csv

Enrichment is only produced when enough significant features are available.

---

## Methods

### MaxQuant filtering
Rows are removed if flagged as:
- Potential contaminant (+)
- Reverse (+)
- Only identified by site (+)

### Feature identifiers
MaxQuant Protein IDs often contains multiple UniProt accessions separated by semicolons.
This pipeline uses the first UniProt accession per row and removes isoform suffixes (for example P12345-2 becomes P12345) as the stable feature identifier.

### Presence filtering
Proteins are filtered to reduce tests dominated by missingness. Default behavior retains proteins meeting a within-condition non-missing threshold (configured in the scripts).

### Normalization and missing-value handling
- Normalization: VSN (vsn::justvsn)
- Imputation: MinProb (MSnbase::impute, method MinProb, default q = 0.01)
- Analysis scale: log2

### Differential abundance testing
- Model: limma
- Design: ~ 0 + condition (one coefficient per condition, no intercept)
- Contrasts: built from metadata/contrasts.csv
- Multiple testing correction: BH (Benjamini-Hochberg)

### Functional enrichment (optional)
GO enrichment is run on significant features after mapping identifiers to gene symbols/Entrez IDs (human mapping is included; adapt as needed for mouse).

---

## Troubleshooting

### Error: mq_column values were not found in proteinGroups
- Open data_raw/proteinGroups.txt and confirm the exact LFQ column headers.
- Copy/paste the exact header strings into metadata/sample_sheet.csv under mq_column.

### Many NAs or non-numeric LFQ values
- Ensure the sample sheet points to LFQ intensity columns, not metadata columns.
- Confirm MaxQuant exported numeric values for those columns.

### PCA clusters by run order / suspected batch effects
- Ensure batch is filled in sample_sheet.csv.
- Extend the limma design to include batch if needed (see R/03_differential_limma.R).

### Enrichment produces no output
- Verify there are enough significant features for the contrast.
- Confirm identifier mapping (UniProt to SYMBOL/ENTREZ) is appropriate for the organism and IDs in your file.

---
