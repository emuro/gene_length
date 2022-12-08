# GeneLength project

## Python3 code (MacBook): 

geneLenght project in the geneLength project dir at enriquem.muro@mac:~/PycharmProjects/geneLength/

----
### get_proteins_ref_proteomes (uniprot, reference proteomes, proteins)
@ get_genes: get_proteins_ref_proteomes.py

It obtains a file for each reference proteome (uniprot). In each file, each entry annotates a 
protein length and the lowest difference in length (aa) with any other protein. The proteins
are sorted by length. 

**main_get_proteins_ref_proteomes.py:** 

* from lib_analysis import constants_analysis as c
* from lib import EM_biopython_extras as bioEx: This parses a Uniprot header

#### Input:
A file with the reference proteomes from Uniprot (19854).
```file_description
For all species,  
@ results/geneLength/outputInputFiles/analysis/some_statistics/stat_description/proteins/  
file: stat_description.protein.uniprot_reference_proteome.tsv
19854 entries + header; superregnum: 
330 archaea, 7997 bacteria, 1588 eukaryota, 9939 viruses = 9915 + 9939 viruses
```

#### Output:
A file with the protein lengths (and uniprot annotations) for each reference proteome obtained from Uniprot.
The proteins are sorted by length (aa).
```file_description
For instance, for homo_sapiens  
@ results/geneLength/ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/  
file: UP000005640_9606.length.tsv
```
----







