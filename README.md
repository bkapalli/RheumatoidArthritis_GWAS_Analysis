# Rheumatoid Arthritis GWAS Analysis  

## Overview  
This repository contains scripts for analyzing Rheumatoid Arthritis (RA) using Genome-Wide Association Studies (GWAS). The analysis is performed using PLINK, EIGENSOFT, METAL, and LD Score Regression

- `quality_control.sh` → Performs initial QC and PCA on genotype data  
- `gwas_analysis.sh` → Runs GWAS using PLINK  
- `meta_analysis.sh` → Conducts meta-analysis with METAL  
- `ldsc_preprocessing.sh` → Prepares summary statistics for LDSC
- `qq_manhattan_plot.R` → Manhattan plots for GWAS results and qq plots for genomic control/meta-analysis results
- `metal.txt` → Configuration file for METAL meta-analysis
