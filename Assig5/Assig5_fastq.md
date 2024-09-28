# Applied Bioinformatics Assignment 5
## Martina Albuja Quintana

### Part 1: Select a genome, then download the corresponding FASTA file.

**Selected genome:** *Halobacterium salinarum* (Strain: 91-R6)

The genome was downloaded and the answers to the questions were obtained using the following script:

https://github.com/Marti-Albuja/MA_AppBioInfo/blob/main/Assig5/Scripts/Assig5_Fasta.sh

Report:

**The size of the file:** 2.3M

**The total size of the genome:** 2429680 bp

**The number of chromosomes in the genome:** 3

**The name (id) and length of each chromosome in the genome:** 
* NZ_CP038631.1 2178608 bp
* NZ_CP038632.1 148406 bp
* NZ_CP038633.1 102666 bp

### Part 2: Generate a simulated FASTQ output for a sequencing instrument of your choice. 

Set the parameters so that your target coverage is 10x.

**Selected genome:** *Halobacterium salinarum* (Strain: 91-R6)

**Simulator:** pbsim3

The genome was downloaded and the answers to the questions were obtained using the following script:

https://github.com/Marti-Albuja/MA_AppBioInfo/blob/main/Assig5/Scripts/Assig5_Fastq_Sim.sh

Report:

**Important:** since pbsim3 generates long read data, the two fastq files generated have different sizes and information

**How many reads have you generated?**

* pbsim_0001.fastq: 2501
* pbsim_0001.fastq: 169

**What is the average read length?**

* pbsim_0001.fastq: 8711 bp
* pbsim_0001.fastq: 8782.6 bp

**How big are the FASTQ files?**

* pbsim_0001.fastq: 42M
* pbsim_0001.fastq: 2.8M

**Compress the files and report how much space that saves**

* pbsim_0001.fastq: 18M
* pbsim_0001.fastq: 1.2M

**Discuss whether you could get the same coverage with different parameter settings (read length vs. read number)**

After trying with different parameters (ex. *--min-length*), I realized that with pbsim3 you have the restriction of not being able to get to a certain coverage by manipulating things like read length and read number. This happens because there is a parameter for *--depth* that controls the coverage you get no matter what you put in the other settings. Even excluding the *--depth* parameter won't work since there is a default depth (20X) that will keep defining the coverage of your simulated reads.

### Part 3: How much data would be generated when covering the Yeast, the Drosophila or the Human genome at 30x?

Now imagine that instead of your genome, each instrument generated reads that cover the Yeast, Drosophila, and Human genomes at 30x coverage (separate runs for each organism). You don't have to run the tool all you need is to estimate.

Using the information you've obtained in the previous points, for each of the organisms, estimate the size of the FASTA file that holds the genome, the number of FASTQ reads needed for 30x, and the size of the FASTQ files before and after compression.

**Note 1:** The fasta file sizes were estimated considering the average genome size of each organism since they are usually really similar.


**Note 2:** If we consider that read lengths are staying the same size as for *Halobacterium salinarum*, then to estimate the other parameters I multiplied the data generated before by 3 (to go from 10X to 30X) and then multiplied by the proportion difference between the size of the genomes.

```
Drosophila (30X): Genome size ~140 Mb 

Size of fasta file: 140M

Number of fastq reads: 

* pbsim_0001.fastq: 435,174
* pbsim_0001.fastq: 29,406

Size of fastq file before compression: 

* pbsim_0001.fastq: 7.3G
* pbsim_0001.fastq: 487M

Size of fastq file after compression:

* pbsim_0001.fastq: 3.1G
* pbsim_0001.fastq: 208.8M

```

```
Yeast (30X): ~12 Mb

Size of fasta file: 12M

Number of fastq reads:

* pbsim_0001.fastq: 37,515
* pbsim_0001.fastq: 2,535

Size of fastq file before compression:

* pbsim_0001.fastq: 630M
* pbsim_0001.fastq: 42M

Size of fastq file after compression:

* pbsim_0001.fastq: 270M
* pbsim_0001.fastq: 18M

```

```
Human (30X): ~2.9 Gb

Size of fasta file: 2.9G

Number of fastq reads:

* pbsim_0001.fastq: 9,063,624
* pbsim_0001.fastq: 612,456 

Size of fastq file before compression:

How big are the FASTQ files?

* pbsim_0001.fastq: 158G
* pbsim_0001.fastq: 10G

Size of fastq file after compression:

* pbsim_0001.fastq: 65G
* pbsim_0001.fastq: 4.3G

```


