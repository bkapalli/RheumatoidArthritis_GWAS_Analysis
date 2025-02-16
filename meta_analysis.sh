module load metal R
# Define input and output directories
# Define PLINK binary file prefix
PLINKFILE=$DATADIR/narac_hg19
metal metal.txt > metal.log
#Extract Necessary Data from the .bim File:
awk '{print $2, $1, $4}' $DATADIR/narac_hg19.bim > $OUTPUTDIR/snp_chr_bp.txt
# Filtering values with significant p values
awk 'NR==1||$10<5e-8{print $0}' merged_meta_data.tbl
# Assuming your output file is named ra_meta_analysis1.tbl
awk '{print $1, $9, $14}' ra_meta_analysis1.tbl > extracted_meta_analysis_data.txt
#Merge CHR and BP Information
awk 'BEGIN {OFS="\t"}
FNR==NR {chr[$1]=$2; pos[$1]=$3; next}
FNR==1 {print $0, "CHR", "BP"}
FNR>1 && $1 in chr {print $0, chr[$1], pos[$1]}
FNR>1 && !($1 in chr) {print $0, "NA", "NA"}' $OUTPUTDIR/snp_chr_bp.txt FS=" " ra_meta_analysis1.tbl FS="\t" >
$OUTPUTDIR/merged_meta_data.tbl
##first, extract chr, bp, and p value from the METAL output:
awk 'BEGIN {FS=OFS="\t"}
NR==1 {print "CHR", "BP", "P-value"} # Print the header for the new file
NR>1 {print $16, $17, $10}' $OUTPUTDIR/merged_meta_data.tbl > $OUTPUTDIR/toplot.txt
Rscript --vanilla metalqqplot.R toplot.txt "Metal QQ Plot"
Rscript --vanilla gwaplot.R toplot.txt "MetaManhattan" metaman.png