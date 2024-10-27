# Applied Bioinformatics Assignment 9
## Martina Albuja Quintana

Use the BAM file generated using reads from the SRA.

In your report, show the commands and the answers for each of these questions:

**Link to Makefile:** 

https://github.com/Marti-Albuja/MA_AppBioInfo/blob/main/Assig9/Makefile

All the following commands can be run in the Makefile by using the following command:

    make analyze

**How many reads did not align with the reference genome?**

To answer this question, the following command was used:

    samtools view -c -f 4 bam/downloaded.bam

The number of reads that did not align to the reference genome was:

3176 reads

**How many primary, secondary, and supplementary alignments are in the BAM file?**

To answer this question, the following commands were used:

    # Primary Alignments
    samtools view -c -F 256 bam/downloaded.bam

    # Secondary Alignments
    samtools view -c -f 256 bam/downloaded.bam

    # Supplementary Alignments
    samtools view -c -f 2048 bam/downloaded.bam

The number of alignments were the following:

**Primary:** 20047 alignments 

**Secondary:** 0 alignments 

**Supplementary:** 47 alignments 

**How many properly-paired alignments on the reverse strand are formed by reads contained in the first pair?**

To answer this question, the following command was used:

    samtools view -c -f 16 -f 2 bam/downloaded.bam

Number of properly-paired alignments on the reverse strand formed by reads contained in the first pair: 8378 alignments 


**Make a new BAM file that contains only the properly paired primary alignments with a mapping quality of over 10**

To create this new bam file using the Makefile, type the command:

    make bam

In the Makefile, the command being used is:

    samtools view -b -q 10 -f 2 -F 256 bam/downloaded.bam > bam/downloaded.filtered.bam

**Compare the flagstats for your original and your filtered BAM file.**

To get the flagstat files using the Makefile, type the command:

    make compare

In the Makefile, the commands being used are:

    samtools flagstat bam/downloaded.bam > bam/downloaded.flagstat
	samtools flagstat bam/downloaded.filtered.bam > bam/downloaded.filtered.flagstat

Original Bam File:

![alt text](https://github.com/Marti-Albuja/MA_AppBioInfo/blob/main/Assig9/Images/Image1.png)

Filtered Bam file:

![alt text](https://github.com/Marti-Albuja/MA_AppBioInfo/blob/main/Assig9/Images/Image2.png)

The % of mapped reads increased from 84.16% to 100% after the filtering process.
