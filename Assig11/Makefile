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
GENOME=Halobacterium.fa

# The name of the GFF file
GFF=Halobacterium.gff

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
	@echo "  analyze: Analyze the alignment"
	@echo "  bam: Create a new BAM file"
	@echo "  compare: Compare the original and filtered BAM files"
	@echo "  variants: Call variants and filter"
	@echo "  vep: Annotate the variants with VEP"

directories:
	mkdir -p ${RDIR} ${PDIR}

reads: directories
	# Download the SRR data
	fastq-dump --split-files -X ${N} -O ${RDIR} ${SRR}

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
	ln -s ncbi_dataset/data/GCF_004799605.1/GCF_004799605.1_ASM479960v1_genomic.fna Halobacterium.fa

	# Make a link to a simpler GFF name
	ln -s ncbi_dataset/data/GCF_004799605.1/genomic.gff Halobacterium.gff

index: 
	# Index the genome for alignment
	bwa index ${GENOME}	

align: 
	# Make the bam directory
	mkdir -p bam

	# Align the reads to the genome
	bwa mem ${GENOME} ${R1} ${R2} | samtools sort -o bam/aligned_reads.bam

	# Index the bam file
	samtools index bam/aligned_reads.bam

stats:  
	# Generate alignment statistics
	samtools flagstat bam/aligned_reads.bam > bam/aligned_reads.flagstat

analyze: 
	# Count the reads that did not align with the reference genome
	@echo "Reads that did not align with the reference genome: $$(samtools view -c -f 4 bam/aligned_reads.bam)"

	# Count number of primary alignments
	@echo "Number of primary alignments: $$(samtools view -c -F 256 -F 2048 bam/aligned_reads.bam)"

	# Count number of secondary alignments
	@echo "Number of secondary alignments: $$(samtools view -c -f 256 bam/aligned_reads.bam)"

	# Count number of supplementary alignments
	@echo "Number of supplementary alignments: $$(samtools view -c -f 2048 bam/aligned_reads.bam)"

	# Count how many properly-paired alignments on the reverse strand are formed by reads contained in the first pair
	@echo "Number of properly-paired alignments on the reverse strand formed by reads contained in the first pair: $$(samtools view -c -f 16 -f 2 bam/aligned_reads.bam)"

bam: 
	# Create new BAM file that contains only properly paired primary alignments with mapping quality of over 10
	samtools view -b -q 10 -f 2 -F 256 bam/aligned_reads.bam > bam/aligned_reads.filtered.bam

compare: 
	# Generate alignment statistics
	samtools flagstat bam/aligned_reads.bam > bam/aligned_reads.flagstat
	samtools flagstat bam/aligned_reads.filtered.bam > bam/aligned_reads.filtered.flagstat

variants:
	# Make the vcf directory
	mkdir -p vcf

	# Call variants
	bcftools mpileup -Ou -f ${GENOME} bam/aligned_reads.filtered.bam | bcftools call -mv -Oz -o vcf/aligned_reads.filtered.vcf.gz

	# Index the vcf file
	bcftools index vcf/aligned_reads.filtered.vcf.gz

	# Count the number of variants
	bcftools stats vcf/aligned_reads.filtered.vcf.gz > vcf/aligned_reads.filtered.vcf.stats

	# Looking for false positives
	bcftools filter -i 'QUAL>=30 && DP>=10 && MQ>=30' vcf/aligned_reads.filtered.vcf.gz > vcf/false_positive_variants.vcf.gz

	bcftools stats vcf/false_positive_variants.vcf.gz > vcf/false_positive_variants.vcf.stats

# VEP is installed in the environment called vep
vep: 
	# Sort and compress the GFF file
	# Needs the double $ to pass the $ from make to bash
	cat ${GFF} | sort -k1,1 -k4,4n -k5,5n -t$$'\t' | bgzip -c > ${GFF}.gz

	# Index the GFF file
	tabix -p gff ${GFF}.gz

	mkdir -p results
	
	micromamba run -n vep \
        ~/src/ensembl-vep/vep \
        -i vcf/aligned_reads.filtered.vcf.gz \
        -o results/vep.txt \
        --gff ${GFF}.gz \
        --fasta ${GENOME} \
        --force_overwrite 

	# Show the resulting files
	ls -lh results/*

all: directories reads download trim index align stats analyze bam compare variants vep

clean:
	rm -rf ${RDIR} ${PDIR} ncbi_dataset.zip Halobacterium.fa* bam vcf results Halobacterium.gff.gz Halobacterium.sorted.gff Halobacterium.gff fastp.html fastp.json README.md md5sum.txt ncbi_dataset Halobacterium.gff.gz.tbi    

.PHONY: all clean directories reads trim download index align stats analyze bam compare variants vep





