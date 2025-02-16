#!/bin/bash
module load R
module load eigensoft
# Load the necessary PLINK module
module load plink/1.90b6.27
# Define input and output directories
DATADIR=
OUTPUTDIR=
# List files in the data directory for verification
echo "Listing all files in the data directory:"
ls -lh $DATADIR
# Define PLINK binary file prefix
PLINKFILE=$DATADIR/narac_hg19
wc -l $DATADIR/narac_hg19.bed
# 754504 /projectnb/bs859/data/RheumatoidArthritis/final_project/narac_hg19.bed
wc -l $DATADIR/narac_hg19.bim
# 544276 /projectnb/bs859/data/RheumatoidArthritis/final_project/narac_hg19.bim
wc -l $DATADIR/narac_hg19.fam
# 2062 /projectnb/bs859/data/RheumatoidArthritis/final_project/narac_hg19.fam
# Check the sex of individuals based on X chromosome genotyping
plink --bfile $PLINKFILE --check-sex --out $OUTPUTDIR/narac_sex_check
# Create a file to exclude the problematic individual
echo "1050200 1050200" > $OUTPUTDIR/exclude.txt
# Step 1: Quality Control with individual exclusion
# Remove SNPs with a high missingness rate (>5%), minor allele frequency (MAF) < 1%, and those not in HWE in
controls (p < 1e-6)
plink --bfile $PLINKFILE --remove $OUTPUTDIR/exclude.txt --maf 0.01 --geno 0.05 --mind 0.05 --hwe 1e-6 --make-
bed --out $OUTPUTDIR/narac_filtered
# Step 3: LD Pruning
# Prune SNPs to reduce the redundancy in the data
plink --bfile $OUTPUTDIR/narac_filtered --indep-pairwise 10000kb 1 0.2 --out $OUTPUTDIR/narac_pruned_data
# Use the pruned data for PCA
plink --bfile $OUTPUTDIR/narac_filtered --extract $OUTPUTDIR/narac_pruned_data.prune.in --make-bed --out
$OUTPUTDIR/narac_pca
#run smartpca
smartpca -p $OUTPUTDIR/smartpca.par > $OUTPUTDIR/smartpca.log
#Plot PCs by case status
# command takes 4 arguments:
# 1) The name of the output file from smartpca
# 2 and 3) the two PCs to plot on x and y axes, respectively
# 4) Number of PCs in the file
# this script assumes the output is from smartpca, so the
# first column is the individual ID and the last column is
# case status (from the plink fam file used to run smartpca)
Rscript --vanilla plotPCs.R narac_pca.evec 1 2 10
Rscript --vanilla plotPCs.R narac_pca.evec 2 4 10
Rscript --vanilla plotPCs.R narac_pca.evec 1 4 10