#Download fasta file from NCBI
datasets download genome accession GCF_004799605.1

# Unzip the data (Overwrite files if they already exist)
unzip -o ncbi_dataset.zip

# Make a link to a simpler name
ln -s ncbi_dataset/data/GCF_004799605.1/GCF_004799605.1_ASM479960v1_genomic.fna Halobacterium.fa

# The location of the genome file
GENOME=Halobacterium.fa

# Sequencing depth
DEPTH=10

# The prefix for the reads
PREFIX=reads/pbsim

# The location of the data error model
MODEL_URL=https://raw.githubusercontent.com/yukiteruono/pbsim3/refs/heads/master/data/QSHMM-RSII.model

# The name of the error model file locally.
MODEL=$(basename ${MODEL_URL})

#
# --- Simulation actions below ---
#

# The read names are created automatically by pbsim3
R1=${PREFIX}_0001.fastq
R2=${PREFIX}_0002.fastq

# Download the data error model
wget -nc ${MODEL_URL}

# Make the reads directory
mkdir -p $(dirname ${PREFIX})

# Run the simulation
pbsim  --strategy wgs --method qshmm --qshmm ${MODEL} \
       --depth ${DEPTH} --genome ${GENOME} --prefix ${PREFIX}

# Statistics on the simulated reads.
seqkit stats ${R1} ${R2}











