# Applied Bioinformatics Assignment 11
## Martina Albuja Quintana

This assignment requires writing a Makefile and a markdown report.

You can reuse the Makefile developed for your previous assignment, which generated a VCF file.

Evaluate the effects of the variants in your VCF file.

Try using a software tool like VEP or snpEff. Add the effect prediction steps to your Makefile and make them part of the workflow.

**Link to Makefile:**

https://github.com/Marti-Albuja/MA_AppBioInfo/blob/main/Assig11/Makefile

**1. Step 1: Download software tool VEP**

```bash
# Create a new conda environment for VEP
micromamba create -y -n vep perl perl-dbi perl-lwp-simple perl-dbd-mysql perl-bio-db-hts

# Activate the environment
micromamba activate vep

# Make a directory for sources
mkdir -p ~/src

# Change to the source directory
cd ~/src

# Clone the VEP repository
git clone https://github.com/Ensembl/ensembl-vep.git

# Change to the VEP directory
cd ensembl-vep

# Install the VEP package
perl INSTALL.pl --NO_HTSLIB --NO_TEST

# Verify the installation
./vep --help

```

**2. Step 2: Modify Makefile with VEP commands**

```bash
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
```

**3. Step 3: Visualize the results**

The results were visualized in the resulting .html file.

In general, 10,181 variants were processed in the genome of *Halobacterium*.

![alt text](https://github.com/Marti-Albuja/MA_AppBioInfo/blob/main/Assig11/Images/Image1.png)

Of all the variants called, 99.3% are SNV or Single Nucleotide Variations/ Polymorphisms and only 0.7% correspond to indels.

![alt text](https://github.com/Marti-Albuja/MA_AppBioInfo/blob/main/Assig11/Images/Image2.png)

Based on the analyzed data, it seems that most of the SNPs cause changes in regions between genes. This could mean that they are probably not changing the transcription or translation of proteins, and therefore, the normal function of the proteins. The other variants are either downstream or upstream of a gene which could potentially be influencing promoter regions.

![alt text](image-18.png)
