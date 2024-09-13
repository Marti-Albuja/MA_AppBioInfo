# Applied Bioinformatics Assigment 2
## Martina Albuja Quintana

### Select an organism and download its corresponding GFF file

#### Greater bamboo lemur (*Prolemur simus*) 

```bash
wget https://ftp.ensembl.org/pub/current_gff3/prolemur_simus/Prolemur_simus.Prosim_1.0.112.gff3.gz
    
gunzip Prolemur_simus.Prosim_1.0.112.gff3.gz
```

### Investigate this file with command line UNIX tools.

### Find answers to the following questions:

* ### Tell us a bit about the organism

The greater bamboo lemur (Prolemur simus) is the most critically endangered lemur species in Madagascar (Wright et al., 2008). Endemic to the island of Madagascar and from the family Lemuridae, it has been listed as one of the 25 most endangered primate species of the world (Hawkins et al., 2018). This is mainly due to deforestation, human hunting, and anthropogenic landscape modification (Wright et al., 2008). Its populations are now fragmented all over Madagascar (Hawkins et al., 2018).
Compared to other lemur species, the greater bamboo lemur is larger with white ear tufts and a broad short muzzle (Hawkins et al., 2018). The greater bamboo lemur is considered a bamboo specialist since 95% of its diet consists of giant or woody bamboo (Wright et al., 2008). Its genome has an approximate size of 2.4 Gb, close to the human's genome (Hawkins et al., 2018). 
  
* ### How many features does the file contain?

#### Typing 
```bash
cat Prolemur_simus.Prosim_1.0.112.gff3 | wc -l
```
#### Prints

```
1387947 features
```

* ### How many sequence regions (chromosomes) does the file contain?

#### Typing
```bash
cat Prolemur_simus.Prosim_1.0.112.gff3 | grep 'sequence-region' | wc -l
```
#### Prints
```
128596 sequence regions
```
* ### How many genes are listed for this organism?

#### Typing
```bash
cat Prolemur_simus.Prosim_1.0.112.gff3 | grep -v '#' > prolemur.gff3
cat prolemur.gff3 | cut -f 3 | sort | uniq -c
```
#### Prints
```
20354 genes
```

* ### What are the top-ten most annotated feature types (col. 3) across the genome?

#### Typing
```bash
cat prolemur.gff3 | cut -f 3 | sort | uniq -c | sort -nr | head
```
#### Prints
```
400695 exon
393577 CDS
128596 region
111882 biological_region
37979 mRNA
20354 gene
5290 ncRNA_gene
1902 transcript
1467 snRNA
 924 rRNA
```

* ### Having analyzed this GFF file, does it seem like a complete and well-annotated organism?

Though the genome is really fragmented (128596 sequence regions), it does seem to be a well-annotated genome. From one side, considering that this animal is a primate species, the number of genes seems to be within the expected number of genes for this family of species (~20,000). Also, exons are usually the most abundant annotated features, so it seems like this genome is following those same patterns.
  
* ### Share any other insights you might note

Considering the endangered status of this animal, this genome and annotation file represents a really important step towards conservation efforts. Though the annotation is pretty good, this genome could still be improved in terms of continuity and assembly. It is interesting to note that this genome has 455 pseudogenes, and it would be an interesting next step if the features "region" and "biological_region" could be better characterized.  

### Github
https://github.com/Marti-Albuja/MA_AppBioInfo

### References
* Hawkins, M. T., Culligan, R. R., Frasier, C. L., Dikow, R. B., Hagenson, R., Lei, R., & Louis, E. E. (2018). Genome sequence and population declines in the critically endangered greater bamboo lemur (*Prolemur simus*) and implications for conservation.BMC genomics,19, 1-15.

* Wright, P. C., Johnson, S. E., Irwin, M. T., Jacobs, R., Schlichting, P., Lehman, S., ... & Tan, C. (2008). The crisis of the critically endangered greater bamboo lemur (*Prolemur simus*).Primate Conservation,23(1), 5-17.

