# HIV Vaccine Microbiome Analysis

## Overview

This repository contains the input data and R scripts used for the analysis presented in:

**Associations between the microbiome and immune responses to an adenovirus-basedBetween the Microbiome and Immune Responses to an Adenovirus-Based HIV-1 candidate vaccine are distinct between Africanafrican and USus cohorts**

Yuhao Li, Daniel J. Stieh, Lindsay Droit, Andrew HyoungJin Kim, Rachel Rodgers, Kathie A. Mihindukulasuriya, Leran Wang, Maria G. Pau, Olive Yua, Herbert W. Virgin, Dan H. Barouch, Megan T. Baldridge, Scott A. Handley

*mSystems* 2026

## Study Description

This study investigates the association between gut bacterial diversity and HIV-1 Ad26-based vaccine immune responses across three geographic regions: the United States (US), Uganda (UG), and Rwanda (RW). 16S rRNA gene amplicon surveys were performed on fecal samples collected from vaccine trial participants at early (week 2) and late (week 26) post-vaccination timepoints.

## Repository Structure

```
├── README.md
├── data/
│   ├── [input phyloseq objects and metadata]
│   └── [output statistics tables]
├── figures/
│   └── [generated figures]
├── scripts/
│   └── [R analysis scripts]
└── [Main_Analysis].Rmd
```

## Input Data

### Phyloseq Objects
The analysis uses phyloseq objects containing:
- 16S rRNA gene amplicon sequence data
- Sample metadata (region, visit, treatment group, etc.)
- Alpha diversity metrics (Richness, Faith's Phylogenetic Diversity)
- Immunological assay data (ELISA, ADCP, ELISpot)

Key data objects:
- `ps1.v9.alpha` - Late timepoint (week 26) with alpha diversity metrics
- `ps1.v2.alpha` - Early timepoint (week 2) with alpha diversity metrics

### Sample Information
- **Regions**: United States (US), Uganda (UG), Rwanda (RW)
- **Treatment groups**: Groups 1-7 (vaccine), Group 8 (placebo)
- **Timepoints**: Week 2 (early), Week 26 (late)

### Immunological Assays
- **ELISA**: Env IgG titers (Clade A, B, C, Con C, Mos1)
- **ADCP**: Antibody-dependent cellular phagocytosis scores
- **ELISpot**: IFNγ responses to Env, Gag, and Pol antigens

## R Scripts and Analysis

### Main Analysis
The primary analysis is contained in the main Rmd file, which includes:
- Data preprocessing and quality control
- Alpha and beta diversity analyses
- Differential abundance testing
- Correlation analyses between bacterial diversity and immune responses
- Figure generation

### Figure 4 Analysis
Figure 4 examines the association between bacterial alpha diversity and vaccine-elicited immune responses using:
- **Panel A**: Spearman's rank-order correlation analysis (summary tables)
- **Panels B, C, D**: ANCOVA linear model comparisons with scatter plots

The Figure 4 code generates:
- Scatter plots with region-specific regression lines
- Multiple statistical output files for review:
  - `Figure4_stats_ANCOVA_pairwise.csv` - Between-region pairwise comparisons
  - `Figure4_stats_Spearman_within_region.csv` - Within-region Spearman correlations
  - `Figure4_stats_Pearson_within_region.csv` - Within-region Pearson correlations
  - `Figure4A_richness_correlations.csv` - Panel A table (Richness)
  - `Figure4A_faithsPD_correlations.csv` - Panel A table (Faith's PD)

## Dependencies

### R Version
R version 4.x or higher recommended

### Required Packages
```r
# Core analysis
library(tidyverse)
library(phyloseq)

# Statistics
library(emmeans)
library(rstatix)

# Visualization
library(ggplot2)
library(cowplot)

# Additional packages as needed
library(vegan)
library(ape)
```

## Running the Analysis

1. Ensure all required packages are installed
2. Set the working directory to the repository root
3. Load the input data (phyloseq objects)
4. Run the main Rmd file to reproduce all analyses and figures

```r
# Example
setwd("/path/to/repository")
rmarkdown::render("Main_Analysis.Rmd")
```

## Output

### Figures
Generated figures are saved to the `figures/` directory in TIFF format at 300 DPI.

### Statistics Tables
Summary statistics and test results are saved to the `data/` directory as CSV files.

## Notes

### Sample Sizes by Region
After excluding placebo (Group 8):
- Rwanda (RW): ~52 samples
- Uganda (UG): ~17 samples  
- United States (US): ~36 samples

The smaller sample size in Uganda may limit statistical power for within-region analyses.

### Figure 4 Statistics
Statistics annotations for Figure 4 panels B, C, D should be reviewed using the output CSV files. Multiple statistical approaches are provided for comparison.

## Contact

Scott A. Handley; shandley@wustl.edu Megan T. Baldridge; mbaldridge@wustl.edu Dan H. Barouch; dbarouch@bidmc.harvard.edu

## License

MIT License

## Citation

If you use this code or data, please cite:

```
[CITATION]
```
