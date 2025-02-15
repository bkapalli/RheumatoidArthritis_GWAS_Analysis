# Rheumatoid Arthritis GWAS Analysis  

## Overview  
This repository contains scripts for analyzing Rheumatoid Arthritis (RA) using Genome-Wide Association Studies (GWAS). The analysis is performed using PLINK, EIGENSOFT, METAL, and LD Score Regression

## Directory Structure  
📂 **scripts/** – Shell scripts for running key analysis steps  
- `quality_control.sh` → Performs initial QC and PCA on genotype data  
- `gwas_analysis.sh` → Runs GWAS using PLINK  
- `meta_analysis.sh` → Conducts meta-analysis with METAL  
- `ldsc_preprocessing.sh` → Prepares summary statistics for LDSC  

📂 **R/** – R scripts for data visualization  
- `plotPCs.R` → Generates PCA plots for population stratification  
- `qq_umich_gc.R` → Creates QQ plots with genomic control  
- `gwaplot.R` → Produces Manhattan plots for GWAS results  
- `metalqqplot.R` → Generates QQ plots for meta-analysis results  

📄 `metal.txt` → Configuration file for METAL meta-analysis  
