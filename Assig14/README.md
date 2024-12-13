# Applied Bioinformatics Assignment 14
## Martina Albuja Quintana

**MakeFile Link:** 

https://github.com/Marti-Albuja/MA_AppBioInfo/blob/main/Assig14/Makefile

To perform a differential expression analysis, please follow the next steps:

We will first be simulating an RNA-Seq library to obtain a design.csv file with the names of each sample and to what group they belong to and a count matrix file (count.csv) with the levels of expression of each gene in each sample. 

**Note:** This procedure uses the *stats* conda environment (make sure to install it before performing the next tasks).

**Step 1: Simulate the RNA-Seq library:**

To do this, use the Makefile target *simulate* 

    make simulate

This target uses the following command:

    src/r/simulate_counts.r -r 5 -n 30000

    # We specify the simulator to create data with 5 replicates and 30,000 genes

It is important to know that the generated data has both upregulated and downregulated genes. The data was simulated using real published data were samples are named as A1, B1, etc and genes as GENE-1, GENE-2, etc.

**Results:**

1. All genes: 30000
2. Genes with data: 7052 
3. Genes that changed: 1500 
4. Changes we can detect: 363 
5. Replicates: 5 
6. Design: design.csv 
7. Counts: counts.csv 

According to the results, there are 363 detectable changes that we will be able to analyze in our dataset.

**Step 2: Identify genes/transcripts that show differential expression:**

To do this, use the Makefile target *identify* 

    make identify

This target uses the following commands:
	
    # Identify differentially expressed genes
    Rscript src/r/edger.r

	# Compare resulting and expected counts of differentially expressed genes 
	Rscript src/r/evaluate_results.r -a counts.csv -b edger.csv

**Results:**

The module edger.r found the following number of statistically significant genes:

1. Significant PVal:  633 ( 9.00 %)
2. Significant FDRs:  259 ( 3.70 %)

633 genes were found to have a p-value lower than 0.5. After applying a correction factor, 256 genes were found to be significant.  

When comparing the counts.csv file and the edger.csv file, we find:

1. 363 in counts.csv 
2. 259 in edger.csv 
3. 248 found in both
4. 115 found only in counts.csv 
5. 11 found only in edger.csv 

Edger was able to identify 248 diferrentially expressed genes that were also marked in the original counts.csv file. Additionally, there were 115 false negatives and 11 false positives gene counts. This means that 115 that were differentially expressed genes were not identified by edger while 11 were misidentified by edger. 

**Step 3: Generate a PCA plot:**

To do this, use the Makefile target *pca* 

    make pca

This target uses the following commands:

    src/r/plot_pca.r -c edger.csv

**Results:**

![alt text](https://github.com/Marti-Albuja/MA_AppBioInfo/blob/main/Assig14/Images/Image1.png)

In the resulting PCA, we can see that the two groups A and B are clearly not correlated with one another and; therefore, are diferentially expressed compared to each other. It is also interesting to note that samples from group B are clustered closer together than samples from group A. These results are explained by 45% of the variance of the data (33% from PCA1 and 12% from PCA2).

**Step 4: Generate a heatmap:**

To do this, use the Makefile target *pca* 

    make heatmake

This target uses the following commands:

    src/r/plot_heatmap.r -c edger.csv

**Results:**

![alt text](image-1.png)

In this image, we have the samples in the X-axis and the genes in the Y-axis. Samples colored red mean that they are upregulated genes while samples colored green mean that they underregulated genes.

Based on all the results obtained, we can clearly see that there are patterns of differential expression between individuals of group A and group B when we analyze the significant genes from our dataset. This is, especially clear, if we observe both the PCA and heatmap.





