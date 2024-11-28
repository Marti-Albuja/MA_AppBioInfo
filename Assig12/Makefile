SHELL := bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

##--------------------------------------------##

# SRR number to download reads
SRR=DRR065898

# Number of reads to download
N=10000

# Downloaded read names
R1=reads/${SRR}_1.fastq
R2=reads/${SRR}_2.fastq

# Trimmed read names
T1=reads/${SRR}_1.trimmed.fastq
T2=reads/${SRR}_2.trimmed.fastq

# The reads directory
RDIR=reads

# The reports directory
PDIR=reports

# The accession number of the genome
ACC=GCA_008636005.1

# The name of the genome file
GENOME=Striga.fa

# The name of the GFF file
GFF=Striga.gff

##--------------------------------------------##

usage:
	@echo "Usage: make <target>"
	@echo ""
	@echo "Targets:"
	@echo "  all: Download reads, trim, align, and analyze"
	@echo "  clean: Remove all generated files"
	@echo "  directories: Create directories for reads and reports"
	@echo "  reads: Download reads and run FastQC"
	@echo "  trim: Trim reads and run FastQC"
	@echo "  download: Download the genome and GFF files"
	@echo "  index: Index the genome for alignment"
	@echo "  align: Align the reads to the genome"
	@echo "  stats: Generate alignment statistics"
	@echo "  variants: Call variants and filter"

directories:
	mkdir -p ${RDIR} ${PDIR}

reads: directories
	# Download the SRR data
	fastq-dump -F -X ${N} -O ${RDIR} --split-files ${SRR}

	# Run fastqc
	fastqc -q -o ${PDIR} ${R1} ${R2}
	
trim: reads
	# Run fastp and trim for quality
	fastp --cut_right -f 30 -F 30 -T 80 -i ${R1} -I ${R2} -o ${T1} -O ${T2}

	# Run fastqc
	fastqc -q -o ${PDIR} ${T1} ${T2}

download:
	# Download fasta file from NCBI
	datasets download genome accession ${ACC} --include genome,gff3

	# Unzip the data (Overwrite files if they already exist)
	unzip -o ncbi_dataset.zip

	# Make a link to a simpler genome name
	ln -sf ncbi_dataset/data/GCA_008636005.1/GCA_008636005.1_SGA_v2.0_genomic.fna Striga.fa

	# Make a link to a simpler GFF name
	ln -sf ncbi_dataset/data/GCA_008636005.1/genomic.gff Striga.gff

index: 
	# Index the genome for alignment
	bwa index ${GENOME}	

align: 
	# Make the bam directory
	mkdir -p bam

	# Align the reads to the genome
	bwa mem ${GENOME} ${R1} ${R2} | samtools sort -o bam/${SRR}.aligned_reads.bam

	# Index the bam file
	samtools index bam/${SRR}.aligned_reads.bam

stats:  
	# Generate alignment statistics
	samtools flagstat bam/${SRR}.aligned_reads.bam > bam/${SRR}.aligned_reads.bam.flagstat

variants:
	# Make the vcf directory
	mkdir -p vcf

	# Call variants
	bcftools mpileup -Ou -f ${GENOME} bam/${SRR}.aligned_reads.bam | bcftools call -mv -Oz -o vcf/${SRR}.aligned_reads.vcf.gz

all: directories reads trim download index align stats variants

clean:
	rm -rf ${RDIR} ${PDIR} ncbi_dataset.zip Striga.fa Striga.gff bam vcf Striga.fa.amb Striga.fa.fai Striga.fa.pac Striga.fa.sa ${GENOME} ${GFF} ${R1} ${R2} ${T1} ${T2} ${SRR}.aligned_reads.bam ${SRR}.aligned_reads.bam.bai ${SRR}.aligned_reads.bam.flagstat vcf/${SRR}.aligned_reads.vcf.gz Striga.fa.ann Striga.fa.bwt README.md fastp.html fastp.json fastp.log fastp.html fastp.json fastp.log md5sum.txt ncbi_dataset

.PHONY: all clean directories reads trim download index align stats variants





