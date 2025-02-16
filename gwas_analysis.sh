#!/bin/bash
module load R
module load plink/1.90b6.27
# Define input and output directories
DATADIR="/projectnb/bs859/data/RheumatoidArthritis/final_project"
OUTPUTDIR="/projectnb/bs859/students/bkapalli/A.project"
# Define PLINK binary file prefix
PLINKFILE=$DATADIR/narac_hg19
##make a temporary file to work on:
awk 'NR>1 {print $0}' $OUTPUTDIR/narac_pca.evec > $OUTPUTDIR/temp1.evec
##remove the colon from the ID:
sed 's/:/\t/g' $OUTPUTDIR/temp1.evec > $OUTPUTDIR/temp2.evec
#add the header:
echo -e "FID IID PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10" > $OUTPUTDIR/narac_pca_covar.txt
awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' $OUTPUTDIR/temp2.evec >>
$OUTPUTDIR/narac_pca_covar.txt
plink --bfile $OUTPUTDIR/narac_filtered --covar $OUTPUTDIR/narac_pca_covar.txt --covar-name PC1-PC10 --out
$OUTPUTDIR/checkPCs --allow-no-sex --logistic no-snp beta
# Calculate allele frequencies for the filtered dataset
plink --bfile $OUTPUTDIR/narac_filtered --freq --out $OUTPUTDIR/allele_freqs
# Merge allele frequencies into the GWAS output
# For female GWAS results
# Perform GWAS for female subjects only, using significant PCs as covariates
plink --bfile $PLINKFILE --filter-females --covar $OUTPUTDIR/narac_pca_covar.txt --covar-name PC2,PC4 --logistic
beta hide-covar --ci 0.95 --out $OUTPUTDIR/GWAS_female
awk 'BEGIN {FS=" "; OFS=" "}
NR==FNR {a[$2]=$3; b[$2]=$4; freq[$2]=$5; next}
FNR==1 {print $0, "ALLELE1", "ALLELE2", "FREQ"}
FNR>1 {if ($2 in a) print $0, a[$2], b[$2], freq[$2]; else print $0, "NA", "NA", "NA"}' $OUTPUTDIR/allele_freqs.frq
$OUTPUTDIR/GWAS_female.assoc.logistic > $OUTPUTDIR/GWAS_female_ready.assoc.logistic
# How many SNPs in this GWAS have p-value < 0.0001
awk 'NR==1||$9<0.0001{print $0}' GWAS_female.assoc.logistic|wc
# Repeat for males
plink --bfile $PLINKFILE --filter-males --covar $OUTPUTDIR/narac_pca_covar.txt --covar-name PC2,PC4 --logistic
beta hide-covar --ci 0.95 --out $OUTPUTDIR/GWAS_male
awk 'BEGIN {FS=OFS=" "}
NR==FNR {allele1[$2]=$3; allele2[$2]=$4; freq[$2]=$5; next}
FNR==1 {print $0, "ALLELE1", "ALLELE2", "FREQ"}
FNR>1 {print $0, allele1[$2], allele2[$2], freq[$2]}' $OUTPUTDIR/allele_freqs.frq
$OUTPUTDIR/GWAS_male.assoc.logistic > $OUTPUTDIR/GWAS_male_ready.assoc.logistic
# How many SNPs in this GWAS have p-value < 0.0001
awk 'NR==1||$9<0.0001{print $0}' GWAS_male.assoc.logistic|wc
# QQ and Manhattan plots using custom R scripts
# QQ plots
Rscript --vanilla $OUTPUTDIR/qq_umich_gc.R $OUTPUTDIR/GWAS_female.assoc.logistic "female GWAS" ADD
Rscript --vanilla $OUTPUTDIR/qq_umich_gc.R $OUTPUTDIR/GWAS_male.assoc.logistic "male GWAS" ADD
# Manhattan plots
Rscript --vanilla $OUTPUTDIR/gwaplot.R $OUTPUTDIR/GWAS_female.assoc.logistic "GWAS Results for Females"
female_manhattan
Rscript --vanilla $OUTPUTDIR/gwaplot.R $OUTPUTDIR/GWAS_male.assoc.logistic "GWAS Results for Males"
male_manhattan
# Summary tables can be created from the PLINK output
awk 'NR==1 || $9<0.05 {print}' $OUTPUTDIR/GWAS_female.assoc.logistic >
$OUTPUTDIR/GWAS_female_summary.txt
awk 'NR==1 || $9<0.05 {print}' $OUTPUTDIR/GWAS_male.assoc.logistic > $OUTPUTDIR/GWAS_male_summary.txt