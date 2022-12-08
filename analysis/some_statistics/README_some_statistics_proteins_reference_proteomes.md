# GeneLength project

## Python3 code (MacBook): 

geneLenght project in the geneLength project dir at enriquem.muro@mac:~/PycharmProjects/geneLength/

----
### some_statistics_proteins_reference_proteomes (genes)
At some_statistics: some_statistics_proteins_reference_proteomes.py  
It gets the statitistical description of the length of the reference proteomes from Uniprot: 
count, mean, std, var, ...  
takes around 1100 seconds ~18.3 minutes

libs,
* lib_analysis:constants_analysis
* pandas
* gzip
* Bio:SeqIO
* numpy
* time:time
* sys
* pathlib

Reference_proteomes from Uniprot. See:
* https://www.uniprot.org/proteomes/?query=*&fil=reference%3Ayes
* ftp://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/reference_proteomes

#### Inputs:
* Reference proteomes,
  
```file_description
path: results/geneLength/outputInputFiles/some_tables/reference_proteomes/reference_proteomes_table_28.5.2021.txt
file: reference_proteomes_table_28.5.2021.txt  
(19854 entries + header; superregnum: 330 archaea, 7997 bacteria, 1588 eukaryota, 9939 viruses  =  
= 9915 + 9939 viruses))
```
* Uniprot annotations (complete reference proteome),  
```file_description
path: data/compressed/ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/ (Bacteria/ Eukaryota/)
file i.e + Eukaryota/UP000005640/UP000005640_9606.fasta.gz
```

#### Output:
The statistical description of the protein length
```file_description
@ results/geneLength/outputInputFiles/analysis/some_statistics/stat_description/proteins/
file: stat_description.protein.reference_proteomes.tsv (also 19854 species)
```

----
**pseudocode:**
1. get the reference proteomes species from *.tsv file

2. For each reference proteome (species):

	* get all the protein sequences annotations: only one protein from each gene
	* using Biopython get the sequences, the sequence length
	* calculate the statistical description of all the protein lengths for the species
	
3. Save the statistical description in the **main output file**  
```file_description
See output
@ results/geneLength/outputInputFiles/analysis/some_statistics/stat_description/proteins/
file: stat_description.protein.uniprot_reference_proteome.tsv
```
----






