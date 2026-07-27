# Proteome and Phosphoproteome Analysis of Rim15 During Cellular Quiescence

This repository contains the R-based analysis pipeline for a quantitative mass spectrometry study examining how the Rim15 kinase shapes the proteome and phosphoproteome of *Saccharomyces cerevisiae* as cells enter quiescence under two different nutrient starvation conditions.

## Background

Rim15 is a conserved serine/threonine kinase that integrates signals from the TORC1, PKA, and Pho80–Pho85 pathways to coordinate entry into quiescence. This project uses SILAC-based quantitative proteomics and phosphoproteomics to compare wild-type (WT) and *rim15Δ* (rim15 knockout) yeast across a time course of:

- **Carbon starvation** (−C)
- **Phosphate starvation** (−Pi)

Four time points (T00, T04, T08, T16) and three biological replicates per condition were collected and analyzed by mass spectrometry. Raw data were processed in MaxQuant before downstream statistical analysis in R.

## Repository Structure

```
.
├── data/
│   ├── proteinGroups.csv          # MaxQuant protein-level output
│   └── Phospho (STY)Sites.csv     # MaxQuant phosphosite-level output
├── tables/
│   ├── master_protein.csv         # Pre-processed protein table (wide format)
│   ├── master_tidy_protein.csv    # Pre-processed protein table (tidy format)
│   ├── normed_ratio_perseus.csv   # Normalized protein ratios (Perseus export)
│   └── normed_pSite_ratio_perseus.csv  # Normalized pSite ratios (Perseus export)
├── plots/                         # Generated PDF figures
├── README_figs/                   # Inline figures for rendered notebooks
│
│  ── Analysis scripts (R Markdown, run in order) ──
├── Data preperation and QC.Rmd        # Step 1 – QC, normalization, PCA, t-tests
├── Proteom analysis.Rmd               # Step 2 – Two-way ANCOVA & KEGG analysis
├── Proteom analysis - ANOVA.Rmd       # Step 3 – One-way ANOVA & WGCNA
├── ANCOVA analysis.Rmd                # Step 4 – Combined proteome/phosphoproteome ANCOVA
├── pSites analysis.Rmd                # Step 5 – Phosphoproteomics analysis
├── Pro-phospho combined.Rmd           # Step 6 – Integrated proteome + phosphoproteome
└── Single protein check.Rmd          # Step 7 – Linearity checks for individual proteins
```

## Analysis Workflow

### Step 1 — Data Preparation and QC (`Data preperation and QC.Rmd`)

- Reads raw `proteinGroups.csv` from MaxQuant.
- Extracts SILAC intensities and ratios, applies the Anders & Huber (2010) normalization for intensities.
- Generates QC boxplots and scatter plots for both raw and normalized data.
- Produces a PCA for each data type.
- Runs t-tests for differential expression at each time point.
- Outputs a cleaned master protein table to `tables/`.

### Step 2 — Proteomics Analysis by ANCOVA (`Proteom analysis.Rmd`)

- Reads the pre-processed master protein table.
- Fits two-way ANCOVA models to test for differential abundance by **genotype** (WT vs. *rim15Δ*) and by **nutrient** condition across the time course.
- Generates volcano plots, heatmaps, and KEGG pathway enrichment visualizations.

### Step 3 — One-Way ANOVA and WGCNA (`Proteom analysis - ANOVA.Rmd`)

- Fits one-way ANOVA to identify proteins significantly changing over time in each condition.
- Runs **WGCNA** (Weighted Gene Co-expression Network Analysis) on normalized ratios for each genotype × nutrient combination to identify co-regulated protein modules.
- Outputs module eigengene PDFs and network topology files (`.RData`).

### Step 4 — Combined ANCOVA (`ANCOVA analysis.Rmd`)

- Integrates the protein group and phosphosite tables.
- Re-runs ANCOVA models and produces summary statistics and visualizations.

### Step 5 — Phosphoproteomics Analysis (`pSites analysis.Rmd`)

- Reads `Phospho (STY)Sites.csv` from MaxQuant.
- Filters sites by localization probability (> 0.7).
- Normalizes phosphosite ratios and runs parallel statistical tests (t-tests, q-values) to identify differentially phosphorylated sites per condition and time point.
- Generates volcano plots, heatmaps, and PCA for phosphosites.

### Step 6 — Integrated Proteome + Phosphoproteome (`Pro-phospho combined.Rmd`)

- Merges log₂ fold-change results from Steps 2 and 5.
- Plots phosphorylation log₂FC against protein abundance log₂FC to identify sites with phosphorylation changes independent of protein level changes (kinase/phosphatase activity).

### Step 7 — Single Protein Linearity Check (`Single protein check.Rmd`)

- Plots normalized ratio and normalized intensity side-by-side for individual proteins (e.g., Igo1) as a sanity check for linearity between the two quantification approaches.

## Dependencies

All analyses are written in **R**. Install the required packages before running the notebooks:

```r
# CRAN packages
install.packages(c(
  "dplyr", "ggplot2", "reshape2", "tidyverse", "GGally",
  "ggpubr", "ggrepel", "gridExtra", "devtools", "ggfortify",
  "limma", "plotly", "broom", "qvalue", "pheatmap",
  "RColorBrewer", "psych", "WGCNA", "gplots", "UpSetR"
))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("org.Sc.sgd.db", "clusterProfiler"))

# GitHub packages
devtools::install_github("omarwagih/rmotifx")
```

## Key Outputs

| File / Directory | Description |
|---|---|
| `plots/ratio_normed_PCA_time.pdf` | PCA of normalized protein ratios colored by time |
| `plots/ratio_normed_heat_all.pdf` | Heatmap of all DE proteins (ratio) |
| `plots/pGroups_volcano.pdf` | Volcano plots for protein groups |
| `plots/pSite_volcano.pdf` | Volcano plots for phosphosites |
| `plots/pSite_ratio_normed_PCA_time.pdf` | PCA of normalized pSite ratios |
| `*_WPM_dynamic(cor>0.7).pdf` | WGCNA module dendrogram (per condition) |
| `ME_*.pdf` | Module eigengene expression profiles |
| `kegg_map_*.pdf` | KEGG pathway maps with DE proteins highlighted |

## Citation / Contact

Analysis by **Siyu Sun**. For questions about the data or methods, please open an issue in this repository.
