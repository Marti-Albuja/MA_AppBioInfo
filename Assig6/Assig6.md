# Applied Bioinformatics Assignment 6
## Martina Albuja Quintana

Write a script that downloads data from the SRALinks to an external site or ENALinks to an external site database and perform a quality control analysis on the data.

When searching for data use the genome you selected for the read simulation assignment as the target organism. Select Illumina or Iontorrent instruments only.

Identify a bad sequencing dataset. You may need to evaluate multiple SRR numbers to find one with poor quality.

**Link to Script:**

https://github.com/Marti-Albuja/MA_AppBioInfo/blob/main/Assig6/Assig6.sh

### Part 1: Quality of downloaded data

**Data Description:** The read data was downloaded from the SRA (NCBI) database for *Halobacterium salinarum* (strain: KBTZ02) pertaining to the BioProject	PRJNA983580. This data has not been published yet as it was submitted on 26-Jun-2023.

**SRA Number:** SRR28572035

This data was run with Illumina MiSeq and produced 4.3M spots/reads corresponding to 2.4G bases.

**Fastqc report for SRR28572035_1.fastq:**

![alt text](image.png)

As seen in the image, reads have really low quality with 3 categories marked as red, 1 in yellow, and 7 in green. In the red categories, we can see that "per base sequence quality" and "per base sequence content" are the ones that seem to have the most problems.

![alt text](image-1.png)

**Fastqc report for SRR28572035_2.fastq:**

![alt text](image-2.png)

In the second set of reads, reads also seem to have really low quality with 2 categories marked as red, 3 in yellow, and 6 in green. In the red categories, we can see that "per base sequence quality" and "per base sequence content" are still the ones that seem to have the most problems.

![alt text](image-3.png)


### Part 2: Quality of filtered data

After filtering the data using the parameters specified in the script (see above), the reads were considerably improved:

**Fastqc report for SRR28572035_1.fastq:**

![alt text](image-4.png)

As we can see in this new image, data was significantly improved. Now there are 9 green categories and 2 yellow categories. Reads are of much better quality and the "per base sequence content improved significantly as well. Additionally, only 826 reads were lost (9174 reads remaining).

![alt text](image-5.png)

**Fastqc report for SRR28572035_2.fastq:**

![alt text](image-6.png)

The same thing happened with the second set of reads. Now there are 9 green categories and 2 yellow categories. Reads are of much better quality and the "per base sequence content" improved significantly as well. Additionally, only 826 reads were lost (9174 reads remaining).

![alt text](image-7.png)


