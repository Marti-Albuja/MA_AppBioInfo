# Applied Bioinformatics Assignment 3
## Martina Albuja Quintana

### Reformat your previous assignment

https://github.com/Marti-Albuja/MA_AppBioInfo/blob/main/Assig2/Assig2.md

### Visualize the GFF file of your choice

### Using a resource of your choice, download the genome and annotation files for an organism of your choice.

#### Organism: *Lupinus albus*

```bash
datasets download genome accession GCA_009771035.1 --include gff3,cds,genome

unzip ncbi_dataset.zip
```

### Use IGV to visualize the annotations relative to the genome

* #### To visualize data in IGV, we first have to index our fasta file

```bash
samtools faidx GCA_009771035.1_CNRS_Lalb_1.0_genomic.fna 
```
#### An overview of annotations in chromosome CM019999.1 (coordinates: CM019999.1:22,108,065-22,109,292) 

**Image 1**: https://github.com/Marti-Albuja/MA_AppBioInfo/blob/main/Assig3/Images/Image1.png

* #### Separate intervals of type "gene" into a different file. If you don't have genes, pick another feature.

```bash
cat genomic.gff | awk ' $3=="gene" { print $0 }' > genes_Lalbus.gff
```
#### An overview of gene annotations vs other features in chromosome CM019999.1 (coordinates: CM019999.1:11,753,971-11,758,884) 

**Image 2**: https://github.com/Marti-Albuja/MA_AppBioInfo/blob/main/Assig3/Images/Image2.png

* #### Using your editor create a GFF that represents intervals in your genome. Load that GFF as a separate track in IGV.

**GFF File Created with 2 Transcripts in Chromosome CM019999.1**: https://github.com/Marti-Albuja/MA_AppBioInfo/blob/main/Assig3/Demo.gff

**Image 3**: https://github.com/Marti-Albuja/MA_AppBioInfo/blob/main/Assig3/Images/Image3.png

