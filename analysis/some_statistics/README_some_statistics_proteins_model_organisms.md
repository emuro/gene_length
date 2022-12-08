# GeneLength project

## Python3 code (MacBook): 

geneLenght project in the geneLength project dir at enriquem.muro@mac:~/PycharmProjects/geneLength/

----
### some_statistics_proteins_model_organisms (genes)
@ some_statistics:some_statistics_proteins_model_organisms.py

libs,
* lib_analysis:constants_analysis
* lib:files
* pandas
* gzip
* Bio:SeqIO
* numpy
* time:time
* sys
* pathlib

It gets the statitistical description of the protein length of some model organisms I have selected manually

#### Inputs:
* Model orgamisms,
```file_description
path: results/geneLength/outputInputFiles/some_tables/model_organisms/proteins/
file: uniprot_well_annotated_organisms.tsv (from P Mier and some that I add like danio rerio)
```

* Uniprot annotations,  
```file_description
path: data/compressed/ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/ (Bacteria/   Eukaryota/)
file i.e + Eukaryota/UP000005640/UP000005640_9606.fasta.gz (Homo Sapiens)
```

#### Output:
The statistical description of the protein length
>@ results/geneLength/outputInputFiles/analysis/some_statistics/stat_description/proteins/  
file: **stat_description.protein.uniprot_model_organism.tsv** (8 entries) 
----
**pseudocode:**
1. get the model organisms from *.tsv file

2. For each model organism (species):

    * Get all the protein sequences annotations, one protein from each gene. For instance, https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz  
    * using Biopython get the sequences, the sequence length
    * calculate the statistical description of the protein lengths for the species
	
3. Save the statistical description in the **main output file**  
```file_description
See output
@ results/geneLength/outputInputFiles/analysis/some_statistics/stat_description/proteins/
   file: stat_description.protein.uniprot_model_organism.tsv
```
----
