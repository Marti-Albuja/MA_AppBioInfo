# Applied Bioinformatics Assignment 12
## Martina Albuja Quintana

To use this Makefile for Automation, please follow the next steps:

**Download the data (SRR) you will be working with:**

    bio search PRJDB4996 -H --csv > design.csv

The species we will be working with is: *Striga asiatica*

**To run the Makefile with the design.csv file, use the following command:**

    # To see the first 10 samples and commands that will be ran

    cat design.csv | head -n 10 | parallel --dry-run --lb -j 4 --colsep , --header : make all SRR={run_accession} SAMPLE={sample_alias}

**Results:**

make all SRR=DRR065898 SAMPLE=SAMD00056127

make all SRR=DRR065899 SAMPLE=SAMD00056127

make all SRR=DRR065902 SAMPLE=SAMD00056127

make all SRR=DRR166235 SAMPLE=SAMD00056127

make all SRR=DRR166236 SAMPLE=SAMD00056127

make all SRR=DRR166237 SAMPLE=SAMD00056127

make all SRR=DRR166240 SAMPLE=SAMD00056127

make all SRR=DRR065894 SAMPLE=SAMD00056127

make all SRR=DRR166239 SAMPLE=SAMD00056127

    # To run the command

    cat design.csv | parallel --lb -j 4 --colsep , --header : make all SRR={run_accession} SAMPLE={sample_alias}

**At the end of the run, you will end up with 21 vcf files:**
DRR065894.aligned_reads.vcf.gz DRR065901.aligned_reads.vcf.gz DRR166239.aligned_reads.vcf.gz
DRR065895.aligned_reads.vcf.gz DRR065902.aligned_reads.vcf.gz DRR166240.aligned_reads.vcf.gz
DRR065896.aligned_reads.vcf.gz DRR065903.aligned_reads.vcf.gz DRR166241.aligned_reads.vcf.gz
DRR065897.aligned_reads.vcf.gz DRR166235.aligned_reads.vcf.gz DRR166242.aligned_reads.vcf.gz
DRR065898.aligned_reads.vcf.gz DRR166236.aligned_reads.vcf.gz DRR166243.aligned_reads.vcf.gz
DRR065899.aligned_reads.vcf.gz DRR166237.aligned_reads.vcf.gz DRR188461.aligned_reads.vcf.gz
DRR065900.aligned_reads.vcf.gz DRR166238.aligned_reads.vcf.gz DRR188462.aligned_reads.vcf.gz

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
**Note 2:** Due to the high amount of computational power needed, it could be possible that the whole pipeline won't be able to run at once. If this happens, run each Makefile target once at a time. 

