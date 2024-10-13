SHELL := bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

usage:
	@echo "Usage: make <target>"
	@echo ""
	@echo "Targets:"
	@echo "  all1: Download and trim reads"
	@echo "  all2: Simulate reads from a reference genome"
	@echo "  directories: Create directories"
	@echo "  download: Download reads"
	@echo "  trim: Trim reads"
	@echo "  genome: Download reference genome"
	@echo "  simulate: Simulate reads"
	@echo "  clean: Remove all files created by targets"

##--------------------------------------------##

# Analysis 1: Download and trim reads 

# SRR number
SRR=SRR28572035

# Number of reads to download
N=10000

# The output read names
R1=reads/${SRR}_1.fastq
R2=reads/${SRR}_2.fastq

# Trimmed read names
T1=reads/${SRR}_1.trimmed.fastq
T2=reads/${SRR}_2.trimmed.fastq

# The reads directory
RDIR=reads

# The reports directory
PDIR=reports


directories:
	mkdir -p ${RDIR} ${PDIR}

download:
	# Download the SRR data
	fastq-dump --split-files -X ${N} -O ${RDIR} ${SRR}

	# Run fastqc
	fastqc -q -o ${PDIR} ${R1} ${R2}

trim:
	# Run fastp and trim for quality
	fastp --cut_right -f 30 -F 30 -T 80 -i ${R1} -I ${R2} -o ${T1} -O ${T2}

	# Run fastqc
	fastqc -q -o ${PDIR} ${T1} ${T2}

all1: directories download trim

##--------------------------------------------##

# Analysis 2: Simulate reads from a reference genome

ACC=GCF_004799605.1

# The location of the genome file
GENOME=Halobacterium.fa

# Sequencing depth
DEPTH=10

# The prefix for the reads
PREFIX=simulate/pbsim

# Read names
RR1=${PREFIX}_0001.fastq
RR2=${PREFIX}_0002.fastq

# The location of the data error model
MODEL_URL=https://raw.githubusercontent.com/yukiteruono/pbsim3/refs/heads/master/data/QSHMM-RSII.model

# The name of the error model file locally.
MODEL=QSHMM-RSII.model

# Download reference genome and simulate reads

genome:
	# Download fasta file from NCBI
	datasets download genome accession ${ACC} 

	# Unzip the data (Overwrite files if they already exist)
	unzip -o ncbi_dataset.zip

	 # Make a link to a simpler name
	ln -s ncbi_dataset/data/GCF_004799605.1/GCF_004799605.1_ASM479960v1_genomic.fna Halobacterium.fa

simulate:
	# Download the data error model
	wget -nc ${MODEL_URL}

	# Make the reads directory
	mkdir -p ${PREFIX}

	# Run the simulation
	pbsim  --strategy wgs --method qshmm --qshmm ${MODEL} --depth ${DEPTH} --genome ${GENOME} --prefix ${PREFIX}

	# Statistics on the simulated reads.
	seqkit stats ${RR1} ${RR2}

all2: genome simulate

##--------------------------------------------##

clean:
	rm -rf ${RDIR} ${PDIR} ${GENOME} ${MODEL} ncbi_dataset.zip Halobacterium.fa QSHMM-RSII.model fastp.json fastp.html ncbi_dataset md5sum.txt simulate 

##--------------------------------------------##
	
.PHONY: clean directories download trim genome simulate all1 all2

all: all1 all2


