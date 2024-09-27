#Download fasta file from NCBI
# Download the genome
datasets download genome accession GCF_004799605.1

# Unzip the data (Overwrite files if they already exist)
unzip -o ncbi_dataset.zip

# Make a link to a simpler name
ln -s ncbi_dataset/data/GCF_004799605.1/GCF_004799605.1_ASM479960v1_genomic.fna Halobacterium.fa

#Extract the file name from simpler name
FILE="Halobacterium.fa"

#
#-----NOTHING NEEDS TO BE CHANGED BELOW THIS LINE-----
#

# Report and count the size of the file
echo "Size of File:$(ls -lh ${FILE} | awk '{print $5}')"

# Report and count the total size of the genome (number of nucleotides)
echo "Genome Size:$(cat ${FILE} | grep -v "^>" | tr -d "\n" | wc -c)"

# Report and count the number of chromosomes in the genome
echo "Number of Chromosomes:$(cat ${FILE} | grep "^>" | wc -l)"

# Report the ID and length of each chromosome in the genome
echo "ID and Length of Each Chromosome:"

#Extract the ID and length of each chromosome in the genome
awk '/^>/ { if (seq) { print id, length(seq) }; id = substr($1, 2); seq = "" } 
     !/^>/ { seq = seq $0 } 
     END { if (seq) { print id, length(seq) } }' ${FILE}













