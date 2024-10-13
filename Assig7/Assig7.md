# Applied Bioinformatics Assignment 7
## Martina Albuja Quintana

In the previous assignments, you were asked to write scripts to simulate reads and in the second assignment you wrote a script to obtain and trim reads for a realistic dataset.

Merge both scripts into a single Makefile.

**Makefile Link:** 

https://github.com/Marti-Albuja/MA_AppBioInfo/blob/main/Assig7/Makefile

**Usage:** This Makefile is divided into two analyses. The analyses can be run separately with:
    
    make all1 # Download and trim reads
    make all2 # Simulate reads from a reference genome

Or they can be run as a whole with:

    make all

In **Analysis 1**, there are three targets to run that will allow you to download Illumina paired-end reads from the SRA NCBI database to trim and improve:

* **directories:** it creates separate directories for downloaded reads and fastqc reports

* **download:** it downloads the reads and runs a prelimimary fastqc analysis to see how reads are before trimming

* **trim:** it trims the paired-end read files using fastp and creates a new fastqc report to see the improvements

In **Analysis 2**, there are two targets to run that will allow you to download a genome from NCBI and simulate long reads using pbsim

* **genome:** downloads and unzips the genome from the NCBI database

* **simulate:** uses the software pbsim and the downloaded genome to simulate long sequencing reads

Additionally, the Makefile contains a target that allows you to clean or eliminate all the files and directories created by any of the two analyses:

    make clean

Finally, all targets in this Makefile have been identified as .PHONY.






