# Applied Bioinformatics Assignment 12
## Martina Albuja Quintana

To use this Makefile for Automation, please follow the next steps:

**Download the data (SRR) you will be working with:**

    bio search PRJNA180 -H --csv > design1.csv

The species we will be working with is: *Desulfurococcus mucosus*

**Due to there being a single-end read file, we will eliminate it from our file using the following command:**

    awk 'NR != 3' design1.csv > design.csv 

**To run the Makefile with the design.csv file, use the following command:**

    # To see the 7 samples and commands that will be ran

    cat design.csv | parallel --dry-run --lb -j 4 --colsep , --header : make all SRR={run_accession} SAMPLE={sample_alias}

**Results:**

make all SRR=SRR064728 SAMPLE=1267

make all SRR=SRR3924292 SAMPLE=4089449

make all SRR=SRR3924293 SAMPLE=4089449

make all SRR=SRR3924294 SAMPLE=4089449

make all SRR=SRR3924295 SAMPLE=4089449

make all SRR=SRR064729 SAMPLE=1267

make all SRR=SRR064727 SAMPLE=1267

    # To run the command

    cat design.csv | parallel --lb -j 4 --colsep , --header : make all SRR={run_accession} SAMPLE={sample_alias}

**At the end of the run, you will end up with 4 vcf files:**

SRR064729.aligned_reads.vcf.gz

SRR3924295.aligned_reads.vcf.gz

SRR3924294.aligned_reads.vcf.gz

SRR3924292.aligned_reads.vcf.gz   

**Note:** There are only 4 output files due to duplication of files in the SRR files.

**Additionally, it is important to know how this Makefile work and what does it do:**

**usage:** to use all the targets of the Makefile use:

    make all

The targets available in the Makefile are:

	clean: Remove all generated files
	directories: Create directories for reads and reports
	reads: Download reads and run FastQC
	trim: Trim reads and run FastQC
	download: Download the genome and GFF files
	index: Index the genome for alignment
	align: Align the reads to the genome
	stats: Generate alignment statistics
	variants: Call variants and filter

**Note:** To download the genome of the species you want to work with replace the ACC code.

    make all ACC=GCF.....

**Once the final .vcf files have been processed, merge and index the files with the following commands:**

    # Merge the vcf files

    bcftools merge SRR064729.aligned_reads.vcf.gz SRR3924292.aligned_reads.vcf.gz SRR3924294.aligned_reads.vcf.gz SRR3924295.aligned_reads.vcf.gz -Oz -o merged.vcf.gz

    # Index the merged file
    bcftools index merged.vcf.gz

**Finally, visualize the different variants in IGV:**

![alt text](https://github.com/Marti-Albuja/MA_AppBioInfo/blob/main/Assig12/Images/Image1.png)

As can be observed in the image, each of the samples has different variants at various positions, but also share some between some of them.
