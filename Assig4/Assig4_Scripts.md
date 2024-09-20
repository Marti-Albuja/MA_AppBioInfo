# Applied Bioinformatics Assignment 4
## Martina Albuja Quintana

### Part 1: Write a Script 

In a previous assignment, you wrote a Markdown report on processing a GFF file to count feature types.

Rewrite all the code from that report as a Bash script. 

**Link to my bash script**

https://github.com/Marti-Albuja/MA_AppBioInfo/blob/main/Assig4/Scripts/GFF_MyScript.sh

Run your script on your original data and verify that it works. 

After running my script, I got the same results as before:

```
Total number of features in the GFF3 file: 1104650

Total number of sequence regions in the GFF3 file: 128596

Top 10 most annotated features in the GFF3 file: 

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

Total number of genes in the GFF3 file: 20354
```

You were also assigned to review someone else's report.

Now, run your script on their data. If the script is reusable, you can replace your variables with theirs and run the script.

I used my script on Robert Bush's data on *Xenopus tropicalis*:

**Link to my bash script when using Robert's data**

https://github.com/Marti-Albuja/MA_AppBioInfo/blob/main/Assig4/Scripts/GFF_Script_Robert.sh


My script was reusable and I was able to obtain results from Robert's data:

```
Total number of features in the GFF3 file: 1492920

Total number of sequence regions in the GFF3 file: 167

Top 10 most annotated features in the GFF3 file: 

598775 exon
557717 CDS
138467 biological_region
68097 five_prime_UTR
52784 three_prime_UTR
49787 mRNA
22107 gene
2023 ncRNA_gene
780 snRNA
535 rRNA

Total number of genes in the GFF3 file: 22107
```

Add more functions to the script that also print some of their results: 

Instead of counting for sequence regions in general, Robert counted the number of chromosomes present in the data which had at least 1000 minimum features. This is why I have adapted the script to count chromosomes with 1000 minimum features instead of sequence regions.

**Link to modified bash script when using Robert's data**

https://github.com/Marti-Albuja/MA_AppBioInfo/blob/main/Assig4/Scripts/GFF_Script_Robert_Updated.sh

Were you able to reproduce their results? Make a note in the report.

I was able to obtain the same 10 chromosomes reported by Robert as well as the same number of features present in each of them: 

```
Total number of chromosomes with at least 1000 features:

195524 1
177280 2
173179 3
169997 4
162333 8
149708 5
138155 6
135136 7
112981 9
77391 10
```
### Part 2: Make use of ontologies

To look for the definitions, parent terms, and children nodes I used:

```
bio explain intron
```

Choose a feature type from the GFF file and look up its definition in the sequence ontology.

**Intron:** A region of a primary transcript that is transcribed, but removed from within the transcript by splicing together the
sequences (exons) on either side of it.

Find both the parent terms and children nodes of the term.

**Parent Terms:** primary_transcript_region 

**Children Nodes**:
* five_prime_intron 
* interior_intron 
* three_prime_intron 
* twintron 
* utr_intron 
* autocatalytically_spliced_intron 
* spliceosomal_intron 
* mobile_intron 
* intron_domain (part_of)
* endonuclease_spliced_intron 
* intronic_regulatory_region (part_of)
* lariat_intron 
* intronic_splicing_silencer (part_of)

Provide a short discussion of what you found.

In general, we can see that the consensus definition of intron refers to the parts of the primary transcript that are cut off before the exons are translated into proteins or RNA sequences. This explains why intron only has primary_transcript_region as its parent term but has a lot of children nodes that are terms related to regions that can be found inside introns. 



