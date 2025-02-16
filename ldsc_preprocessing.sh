#!/bin/bash
# Define directories and files
RA_DIR='/projectnb/bs859/students/bkapalli/A.project/Okada_2014'
LDSCORES_DIR='/projectnb/bs859/data/ldscore_files'
OUTPUT_DIR='/projectnb/bs859/students/bkapalli/A.project'
module load R
module load python2
module load ldsc
# Preprocessing GWAS summary statistics with log transformation
zcat $RA_DIR/RA_GWASmeta_Asian_v2.txt.gz | \
awk 'BEGIN{OFS="\t"} {if(NR==1) print $0, "log_OR_A1"; else print $0, log($6)}' | \
gzip > $RA_DIR/RA_GWASmeta_Asian_v2_with_log.gz
zcat $RA_DIR/RA_GWASmeta_European_v2.txt.gz | \
awk 'BEGIN{OFS="\t"} {if(NR==1) print $0, "log_OR_A1"; else print $0, log($6)}' | \
gzip > $RA_DIR/RA_GWASmeta_European_v2_with_log.gz
# Prepare summary statistics by formatting them correctly with munge_sumstats.py
# European ancestry
# Total sample size = 14361 (cases) + 43923 (controls)
munge_sumstats.py \
--sumstats $RA_DIR/RA_GWASmeta_European_v2_with_log.gz \
--snp SNPID \
--a1 A1 \
--a2 A2 \
--signed-sumstats log_OR_A1,0 \
--p P-val \
--N 58284 \
--merge-alleles $LDSCORES_DIR/w_hm3.snplist \
--out $OUTPUT_DIR/RA_EUR
# Asian ancestry
# Total sample size = 4873 (cases) + 17642 (controls)
munge_sumstats.py \
--sumstats $RA_DIR/RA_GWASmeta_Asian_v2_with_log.gz \
--snp SNPID \
--a1 A1 \
--a2 A2 \
--signed-sumstats log_OR_A1,0 \
--p P-val \
--N 22515 \
--merge-alleles $LDSCORES_DIR/w_hm3.snplist \
--out $OUTPUT_DIR/RA_ASN
# Estimate heritability using LD score regression
# European ancestry
ldsc.py \
--h2 $OUTPUT_DIR/RA_EUR.sumstats.gz \
--ref-ld $LDSCORES_DIR/UKBB.ALL.ldscore/UKBB.EUR.rsid \
--w-ld $LDSCORES_DIR/UKBB.ALL.ldscore/UKBB.EUR.rsid \
--out $OUTPUT_DIR/RA_h2_EUR
# East Asian
ldsc.py \
--h2 $OUTPUT_DIR/RA_ASN.sumstats.gz \
--ref-ld $LDSCORES_DIR/UKBB.ALL.ldscore/UKBB.EAS.rsid \
--w-ld $LDSCORES_DIR/UKBB.ALL.ldscore/UKBB.EAS.rsid \
--out $OUTPUT_DIR/RA_h2_ASN_EAS
# Central/South Asian
ldsc.py \
--h2 $OUTPUT_DIR/RA_ASN.sumstats.gz \
--ref-ld $LDSCORES_DIR/UKBB.ALL.ldscore/UKBB.CSA.rsid \
--w-ld $LDSCORES_DIR/UKBB.ALL.ldscore/UKBB.CSA.rsid \
--out $OUTPUT_DIR/RA_h2_ASN_CSA
# Middle Eastern
ldsc.py \
--h2 $OUTPUT_DIR/RA_ASN.sumstats.gz \
--ref-ld $LDSCORES_DIR/UKBB.ALL.ldscore/UKBB.MID.rsid \
--w-ld $LDSCORES_DIR/UKBB.ALL.ldscore/UKBB.MID.rsid \
--out $OUTPUT_DIR/RA_h2_ASN_MID
done
echo "All analyses are complete."