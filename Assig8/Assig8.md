# Applied Bioinformatics Assignment 8
## Martina Albuja Quintana

Use the Makefile developed for your previous assignment that contains rules for simulating reads and obtaining reads from SRA. Make a new version of the Makefile that includes rules for aligning reads to a reference genome.

**Makefile Link:** 

https://github.com/Marti-Albuja/MA_AppBioInfo/blob/main/Assig8/Makefile

**Note:** Compared to the Makefile created for the previous assigment, I have improved it so that it is no longer separated in sections but is easier to read as a whole. 

In this way, to run all the targets you would only write:

    make all

**Add new targets called index and align to the Makefile.**

The index target should create an index for the reference genome. The align target should align the reads to the reference genome. The output of the align target should be a sorted, indexed BAM alignment file.

Apart from the targets that were already included in the previous Makefile:

* **directories:** it creates separate directories for downloaded reads and fastqc reports

* **download:** it downloads the reads and runs a prelimimary fastqc analysis to see how reads are before trimming

* **trim:** it trims the paired-end read files using fastp and creates a new fastqc report to see the improvements

* **genome:** downloads and unzips the genome from the NCBI database

* **simulate:** uses the software pbsim and the downloaded genome to simulate long sequencing reads

This new Makefile includes three new targets:

* **index:** indexes the reference genome to be able to visualize it in IGV.

* **align:** aligns both the downloaded and simulated reads to the genome, sorts and indexes the resulting bam file to visualize in IGV.

* **stats:** generates alignment statistics for both downloaded and simulated bam files.

**Note:** As the read simulator I used (pbsim), did not create paired reads, the script was modified to map the two sets of simulated reads separately, index them, and then merge them into a single file for better visualization. 


**Visualize the resulting BAM files for your simulated reads and the reads downloaded from SRA.**

![alt text](https://github.com/Marti-Albuja/MA_AppBioInfo/blob/main/Assig8/Images/Image1.png)

**Generate alignment statistics for the reads from both sources, simulated and SRA.**

Downloaded reads alignment statistics:

![alt text](https://github.com/Marti-Albuja/MA_AppBioInfo/blob/main/Assig8/Images/Image2.png)

Simulated reads alignment statistics:

![alt text](image-10.png)

**Briefly describe the differences between the two datasets.**

As can be seen in the visualization with IGV and the % of reads mapped to the genome, the simulated reads were the ones that best aligned and covered the regions of the reference genome for *Halobacterium salinarum*. It is important to note that the downloaded reads are short paired illumina reads while the simulated reads are long single paired reads. 

In the downloaded reads, we can see that reads do map to the genome but are scattered and do not cluster together as much, giving a really low coverage. This can also be seen in the statistics were the % amount of reads that map to the genome is 84.16%.

In contrast, the simulated reads have a 96.66% of mapped reads to the reference genome with a higher coverage observed in the IGV graph. Reads are colored based on their direction (forward or reverse) but we cannot pair them since they are single-end reads. Nonetheless, they align much better than the downloaded reads. 













