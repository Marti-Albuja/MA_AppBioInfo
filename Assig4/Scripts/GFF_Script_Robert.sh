#Set the URL to the file you want to download
URL="https://ftp.ensembl.org/pub/current_gff3/xenopus_tropicalis/Xenopus_tropicalis.UCB_Xtro_10.0.112.gff3.gz"

#Extract the file name from the URL 
FILE=$(basename ${URL})

#Unzipped file name
UNZIPPED=$(basename ${URL} .gz) 

#File name for the extracted gene features
GENES="genes.gff3"

#
#-----NOTHING NEEDS TO BE CHANGED BELOW THIS LINE-----
#

#Download the GFF3 file
wget ${URL}

#Unzip the file
gunzip ${FILE}

# Report and count how many features are in the GFF3 file
echo "Total number of features in the GFF3 file: $(cat ${UNZIPPED} | grep -v "^#" | wc -l)"

#Report and count how many sequence regions are in the GFF3 file
echo "Total number of sequence regions in the GFF3 file: $(cat ${UNZIPPED} | grep 'sequence-region' | wc -l)"

#Print the top 10 most annotated features in the GFF3 file
echo "Top 10 most annotated features in the GFF3 file: $(cat ${UNZIPPED} | grep -v "^#" | awk '{print $3}' | sort | uniq -c | sort -nr | head -n 10)"

#Extract the gene features from the GFF3 file
cat ${UNZIPPED} | grep -v "^#" | awk '$3 == "gene"' > ${GENES}

#Report and count how many genes are in the GFF3 file
echo "Total number of genes in the GFF3 file: $(cat ${GENES} | wc -l)"













