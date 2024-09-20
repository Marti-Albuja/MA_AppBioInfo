#Set the URL to the file you want to download
URL="https://ftp.ensembl.org/pub/current_gff3/prolemur_simus/Prolemur_simus.Prosim_1.0.112.gff3.gz"

#Extract the file name from the URL 
FILE=$(basename ${URL})

#File name for the extracted gene features
GENES="genes.gff3"

#
#-----NOTHING NEEDS TO BE CHANGED BELOW THIS LINE-----
#

#Download the GFF3 file
wget ${URL}

#Unzip the file
gunzip ${FILE}

#Count how many features are in the GFF3 file
cat ${FILE} | grep -v "^#" ${FILE} | wc -l

#Count how many sequence regions are in the GFF3 file
cat ${FILE} | grep 'sequence-region' | wc -l

#Print the top 10 most annotated features in the GFF3 file
cat ${FILE} | grep -v "^#" | awk '{print $3}' | sort | uniq -c | sort -nr | head -n 10

#Extract the gene features from the GFF3 file
cat ${FILE} | grep -v "^#" | awk '$3 == "gene"' > ${GENES}

#Count how many genes are in the GFF3 file
cat ${GENES} | wc -l











