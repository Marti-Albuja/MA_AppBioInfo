# Applied Bioinformatics Assignment 13
## Martina Albuja Quintana

**MakeFile Link:** 

https://github.com/Marti-Albuja/MA_AppBioInfo/blob/main/Assig13/Makefile

To generate an RNA-Seq count matrix, please follow the next steps:

We will be working with RNA-Seq data from *Lucilia cuprina*, with five files of data generated with Illumina NovaSeq 6000 to study how larval and adult stages of *Lucilia cuprina* employ different sensory systems.

**Download the RNA-Seq data (SRR) you will be working with:**

We will be using the following csv file: https://github.com/Marti-Albuja/MA_AppBioInfo/blob/main/Assig13/sample.csv

The data will be downloaded using the target *reads* from the Makefile using a parallel command:

    cat sample.csv | parallel --lb -j 4 --colsep , --header : make reads SRR={SRR}

**Download the genome and gff3 file for Lucila cuprina:**

To do this, use the Makefile target *download*

    make download

**Before running HiSat2 to align RNA-Seq data to a reference genome, we will have to index the reference genome:**

To do this, use the Makefile target *index*

    make index

This target uses the follwoing command:

    hisat2-build ${GENOME} ${GENOME}-index

**We will now run HiSat2 to allign the RNA-Seq data to the reference genome and create a sorted, indexed bam file for each sample used:**

To do this, use the Makefile target *align*

    make align

To do this for all the data sets we have, use the parallel command:

    cat sample.csv | parallel --lb -j 4 --colsep , --header : make align SRR={SRR}

**Note:** In this step, we will also create a bigwig coverage file and its index using the command bamCoverage from the deeptools package.

**Results:**

We will now visualize the alignment of the 5 RNA-Seq files in IGV:

![alt text](image-2.png)

From this example of a part of the genome, we can say that in the 5 files used there are sites where the expression levels are high in every single file and sites that do no have expression at all. Then, there are also sites that are expressed in one or only a few of the RNA-Seq Files.

**Finally, we will create an RNA-Seq count matrix for the five files:**

To do this, we will use the following command:

    featureCounts -a Lucila.gff -o matrix.txt SRR31641679_aligned_reads.bam SRR31641680_aligned_reads.bam SRR31641681_aligned_reads.bam SRR31641682_aligned_reads.bam SRR31641683_aligned_reads.bam 

From the summary.txt matrix file, we see that each of the five different RNA-Seq files has simiilar number of reads alligned to the reference genome; and therefore, probably, similar expression levels.

SRR31641679_aligned_reads.bam = 13987
SRR31641680_aligned_reads.bam =	14172
SRR31641681_aligned_reads.bam = 14970
SRR31641682_aligned_reads.bam = 14852
SRR31641683_aligned_reads.bam = 14851

**Note:** looking at read counts of individuals regions was very difficult because the matrix.txt file created had an enormous amount of information.




