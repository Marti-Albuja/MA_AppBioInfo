SHELL := bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

##--------------------------------------------##

# SRR number to download reads
SRR=SRR28572035

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

ACC=GCF_004799605.1

# The name of the genome file
GENOME= Halobacterium.fa

# Sequencing depth
DEPTH=10

# Prefix for the simulated reads
PREFIX=simulate/pbsim

# Simulated Read Names
RR1=${PREFIX}_0001.fastq
RR2=${PREFIX}_0002.fastq

# The location of the data error model
MODEL_URL=https://raw.githubusercontent.com/yukiteruono/pbsim3/refs/heads/master/data/QSHMM-RSII.model

# The name of the error model file locally.
MODEL=QSHMM-RSII.model

##--------------------------------------------##

usage:
	@echo "make all: Download, trim, align, and analyze the data"
	@echo "make clean: Remove all downloaded files"
	@echo "make directories: Create the directories for the data"
	@echo "make download: Download the data"
	@echo "make trim: Trim the data"
	@echo "make genome: Download the genome"
	@echo "make simulate: Simulate reads"
	@echo "make index: Index the genome"
	@echo "make align: Align the reads to the genome"
	@echo "make stats: Generate alignment statistics"
	@echo "make analyze: Analyze the alignment"
	@echo "make bam: Create a new BAM file"
	@echo "make compare: Compare the alignment statistics"
	@echo "make variants: Call variants"

directories:
	mkdir -p ${RDIR} ${PDIR}

download: directories
	# Download the SRR data
	fastq-dump --split-files -X ${N} -O ${RDIR} ${SRR}

	# Run fastqc
	fastqc -q -o ${PDIR} ${R1} ${R2}

trim: download
	# Run fastp and trim for quality
	fastp --cut_right -f 30 -F 30 -T 80 -i ${R1} -I ${R2} -o ${T1} -O ${T2}

	# Run fastqc
	fastqc -q -o ${PDIR} ${T1} ${T2}

genome: 
	# Download fasta file from NCBI
	datasets download genome accession ${ACC} 

	# Unzip the data (Overwrite files if they already exist)
	unzip -o ncbi_dataset.zip

	# Make a link to a simpler name
	ln -s ncbi_dataset/data/GCF_004799605.1/GCF_004799605.1_ASM479960v1_genomic.fna Halobacterium.fa

simulate: genome
	# Download the data error model
	wget -nc ${MODEL_URL}

	# Make the reads directory
	mkdir -p ${PREFIX}

	# Run the simulation
	pbsim  --strategy wgs --method qshmm --qshmm ${MODEL} --depth ${DEPTH} --genome ${GENOME} --prefix ${PREFIX}

	# Statistics on the simulated reads.
	seqkit stats ${RR1} ${RR2}

index: genome
	# Index the genome for alignment
	bwa index ${GENOME}	

align: trim genome simulate
	# Make the bam directory
	mkdir -p bam

	# Align the downloaded reads to the genome
	bwa mem ${GENOME} ${R1} ${R2} | samtools sort -o bam/downloaded.bam

	# Index the bam file
	samtools index bam/downloaded.bam

	# Align the first set of simulated reads to the genome
	bwa mem ${GENOME} ${RR1} | samtools sort -o bam/simulated_1.bam

	# Index the bam file for the first set of simulated reads
	samtools index bam/simulated_1.bam

	# Align the second set of simulated reads to the genome
	bwa mem ${GENOME} ${RR2} | samtools sort -o bam/simulated_2.bam

	# Index the bam file for the second set of simulated reads
	samtools index bam/simulated_2.bam

	# Merge the two sets of simulated reads
	samtools merge bam/simulated.bam bam/simulated_1.bam bam/simulated_2.bam

	# Index the merged bam file
	samtools index bam/simulated.bam

stats: 
	# Generate alignment statistics
	samtools flagstat bam/downloaded.bam > bam/downloaded.flagstat
	samtools flagstat bam/simulated.bam > bam/simulated.flagstat

analyze: 
	# Count the reads that did not align with the reference genome
	@echo "Reads that did not align with the reference genome: $$(samtools view -c -f 4 bam/downloaded.bam)"

	# Count number of primary alignments
	@echo "Number of primary alignments: $$(samtools view -c -F 256 bam/downloaded.bam)"

	# Count number of secondary alignments
	@echo "Number of secondary alignments: $$(samtools view -c -f 256 bam/downloaded.bam)"

	# Count number of supplementary alignments
	@echo "Number of supplementary alignments: $$(samtools view -c -f 2048 bam/downloaded.bam)"

	# Count how many properly-paired alignments on the reverse strand are formed by reads contained in the first pair
	@echo "Number of properly-paired alignments on the reverse strand formed by reads contained in the first pair: $$(samtools view -c -f 16 -f 2 bam/downloaded.bam)"

bam: 
	# Create new BAM file that contains only properly paired primary alignments with mapping quality of over 10
	samtools view -b -q 10 -f 2 -F 256 bam/downloaded.bam > bam/downloaded.filtered.bam

compare: 
	# Generate alignment statistics
	samtools flagstat bam/downloaded.bam > bam/downloaded.flagstat
	samtools flagstat bam/downloaded.filtered.bam > bam/downloaded.filtered.flagstat

variants: 
	# Call variants
	bcftools mpileup -Ou -f ${GENOME} bam/downloaded.filtered.bam | bcftools call -mv -Ob -o bam/downloaded.filtered.vcf.gz

	# Index the vcf file
	bcftools index bam/downloaded.filtered.vcf.gz	

	# Count the number of variants
	bcftools stats bam/downloaded.filtered.vcf.gz > bam/downloaded.filtered.vcf.stats

	# Looking for false positives
	bcftools filter -i 'QUAL>=30 && DP>=10 && MQ>=30' bam/downloaded.filtered.vcf.gz > false_positive_variants.vcf.gz

	bcftools stats false_positive_variants.vcf.gz > false_positive_variants.vcf.stats

all: directories download trim genome simulate index align stats analyze bam compare variants

clean:
	rm -rf ${RDIR} ${PDIR} ${GENOME} ${MODEL} ncbi_dataset.zip Halobacterium.fa QSHMM-RSII.model fastp.json fastp.html ncbi_dataset md5sum.txt simulate README.md ${PREFIX} ${GENOME}.bwt ${GENOME}.pac ${GENOME}.ann ${GENOME}.sa ${GENOME}.amb ${GENOME}.fai ${GENOME}.dict bam   

.PHONY: clean directories download trim genome simulate all usage index align stats analyze bam compare variants





