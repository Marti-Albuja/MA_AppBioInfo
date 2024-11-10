# Applied Bioinformatics Assignment 10
## Martina Albuja Quintana

You can reuse the Makefile developed for your previous assignment, which generated a BAM file from SRA reads. You may need to get more data to obtain sufficient coverage over the genome. If the data shows no variants, find another dataset that does.

Call variants on the BAM file and discuss some information on the variants you found.

Your eye is an excellent variant caller.

**Link to Makefile:** 

https://github.com/Marti-Albuja/MA_AppBioInfo/blob/main/Assig10/Makefile

To obtain the variants the following commands were added to the Makefile:

    variants: 
    
    # Call variants
	bcftools mpileup -Ou -f ${GENOME} bam/downloaded.filtered.bam | bcftools call -mv -Ob -o bam/downloaded.filtered.vcf.gz

	# Index the vcf file
	bcftools index bam/downloaded.filtered.vcf.gz	


Some stats were also calculated for the VCF file: 

    bcftools stats bam/downloaded.filtered.vcf.gz > bam/downloaded.filtered.vcf.stats

**Number of variants found:** 10,181 variants

To analyze some of the variants called, the following command was used:

     bcftools view downloaded.filtered.vcf.gz | less

**Examples of variants found:**

**Example 1:**

Heading of VCF File:

    #CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  


SNP:

    NZ_CP038631.1  511  .  A  G  8.99921 .  DP=1;SGB=-0.379885;MQ0F=0;AC=2;AN=2;DP4=0,0,0,1;MQ=60  GT:PL   1/1:38,3,0

SNP:

    NZ_CP038631.1   568  .  A  G  8.99921 . DP=1;SGB=-0.379885;MQ0F=0;AC=2;AN=2;DP4=0,0,1,0;MQ=60   GT:PL  1/1:38,3,0

SNP: 

    NZ_CP038631.1   795  .  C  T  68.4148 .  DP=3;VDB=0.28;SGB=-0.453602;MQSBZ=0;MQ0F=0;AC=2;AN=2;DP4=0,0,1,1;MQ=60 GT:PL   1/1:98,6,0

1. In general, there is a great number of SNP variants. As seen in the examples above, we can see nucleotide specific changes from A to G and from C to T.
2. We can also see that the quality score is similar in the first two entries, but much higher (more confidence) in the third entry. 
3. The Depth Coverage of the two first variants is 1 while for the third is 3.

**Example 2:**

Indel:

    NZ_CP038631.1 2138 . AGGGGGGG AGGGGGG 31.4175 . INDEL;IDV=6;IMF=1;DP=6;VDB=0.922588;SGB=-0.590765;MQSBZ=0; MQ0F=0;AC=2;AN=2;DP4=0,0,2,3;MQ=60 GT:PL 1/1:61,15,0

1. There are also indels in the variants called from this data set. Here, we see a deletion of a G in the analyzed sequence compared to the reference.
2. The quality of this variant is relatively high (31.4175).
3. The Depth Coverage (6) is higher than the ones observed in the SNP variants analyzed above.

**Example 3:**

Indel:

    NZ_CP038631.1   22417  .  CTG  C   5.04449 .  INDEL;IDV=2;IMF=1;DP=2;VDB=0.92;SGB=-0.453602;MQSBZ=0; MQ0F=0;AC=2;AN=2;DP4=0,0,1,1;MQ=60  GT:PL 1/1:33,6,0

1. This is another example of a deletion in position 22,417.
2. In this case, the quality score (5.04449) is lower than all the other variants analyzed previously.
3. The Depth Coverage is 2 which means there are 2 reads that support this variant.

**Verify the variant caller's results by looking at a few example alignments in the BAM file.**

![alt text](https://github.com/Marti-Albuja/MA_AppBioInfo/blob/main/Assig10/Images/Image1.png)

As we can see in the image, there are a lot of variants (especially SNPs) that can be identified in this particular region of the genome in different shades of blue. 

**Find examples where the variant caller did not work as expected: false positives, false negatives, etc.**

**False Positives:** A false positive occurs when the variant caller detects a variant that doesn’t exist.

To look for False Positives the following commands were used:

    bcftools filter -i 'QUAL>=30 && DP>=10 && MQ>=30' bam/downloaded.filtered.vcf.gz > false_positive_variants.vcf.gz

	bcftools stats false_positive_variants.vcf.gz > false_positive_variants.vcf.stats

False Positive variants were considered variants with a quality score below 30, a depth coverage below 10, and an MQ below 30.

**False Positives Identified:** 23

**False Negatives:** A false negative occurs when a true variant is missed during the variant calling process.

To look for False Negatives, we assume that these may occur in areas with low coverage and that is why they were not detected. In this way, I used a more visual approach and scanned the alignments in IGV to look for potential regions where False Negatives could occur.

**Example 1:**

![alt text](https://github.com/Marti-Albuja/MA_AppBioInfo/blob/main/Assig10/Images/Image2.png)

**Example 2:**

![alt text](https://github.com/Marti-Albuja/MA_AppBioInfo/blob/main/Assig10/Images/Image3.png)


