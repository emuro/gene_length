# gene_length
## The emergence of eukaryotes signaled by a critical protein-to-gene length ratio

This repository contains the data and programs needed to reproduce the results reported 
in our article; we start explaining how to obtain the data necessary to produce those 
results.  

**The structure of this repository:**  
 - **README.md** This file, that guides you all over the repository
 - **main_tables** where the tables needed to reproduce our results is 
 - **main_work** contains the software needed to reproduce the main results/figures of our work is. This directory contains **main_suppl**, where the programs for the supplementary material are.

---
### Data: the annotations were downloaded from public repositories:

#### Genes
The protein coding gene annotations were obtained from different Ensembl's webservers 
for [prokaryotes (archaea, bacteria)](https://bacteria.ensembl.org), 
[protists](https://protists.ensembl.org), [plants](https://plants.ensembl.org), 
[fungi](https://fungi.ensembl.org), [metazoa](https://metazoa.ensembl.org), 
[vertebrates](https://www.ensembl.org).  


| Ensembl ftp site by Kingdom/division                                          | Release            |  
| :---------------------------------------------------------------------------  | :----------------- |  
| [prokaryotes: archaea, bacteria](http://ftp.ensemblgenomes.org/pub/bacteria/) | ensemblgenomes 49  |  
| [protists](http://ftp.ensemblgenomes.org/pub/protists/)                       | ensemblgenomes 49  |  
| [plants](http://ftp.ensemblgenomes.org/pub/plants/)                           | ensemblgenomes 49  |  
| [fungi](http://ftp.ensemblgenomes.org/pub/fungi/)                             | ensemblgenomes 49  |  
| [metazoa](http://ftp.ensemblgenomes.org/pub/metazoa/)                         | ensemblgenomes 49  |  
| [vertebrates](https://ftp.ensembl.org/pub/)                                   | ensembl 98         |  

The gzip compressed *.gtf.gz (General Transfer Format) gene annotation files were downloaded 
for the different species preserving the directories' structure of the FTP Ensembl 
repositories. For instance, the annotation file of _Homo sapiens_ was downloaded at 
```
Homo_sapiens.GRCh38.98.gtf.gz @
our_mnt_dir + data/compressed/ + "ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/"
```
our_mnt_dir is the local directory where the data was downloaded.


#### Proteins
[Reference proteomes](https://www.uniprot.org/proteomes/?query=*&fil=reference%3Ayes)
were downloaded from [Uniprot](https://www.uniprot.org/), 
each has a unique Uniprot-identifier (UPID). 
A [description](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/README) 
of the reference proteomes, as well as a [table](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/README) 
associating UPIDs, taxonomy_ids, species names, etc is available at Uniprot.

The reference proteomes for the different taxonomical divisions provided by Ensembl (Viruses were not considered) were downloaded from 
[Uniprot FTP repository](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/) on 28.5.2021. 
For each species, a fasta file containing its reference proteome was downloaded. 
The directory structure of the FTP repository was preserved.   
For instance, for _Homo sapiens_ (UPID: UP000005640 and taxonomy id:9606): 
```
UP000005640_9606.fasta.gz @
our_mnt_dir + /data/compressed/ + "ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/"
```


---
### main_tables
For protein coding genes, proteins, and the intersection set between them (merged). The files (*.tsv) are provided in standard [tab-separated values](https://en.wikipedia.org/wiki/Tab-separated_values).
- stat_protCodGenes.tsv (one header line + 33,629 entries)
- stat_proteins.tsv (one header line + 9,915 entries)
- stat_merged (one header line + 6,521 entries)


#### The files contains the next number of entries per taxonomical division:
stat_protCodGenes.tsv (one header line + 33,629 entries):
| counts | regnum               |  
|-----:  |:----------           |
| 31943  | bacteria<sup>*</sup> |
| 237    | protists    |
| 96     | plants      |
| 1014   | fungi       |
| 115    | metazoa     |
| 224    | vertebrates |
33629 entries in total  

<sup>*</sup>In the annotation from Ensembl bacteria includes also archaea.

stat_proteins.tsv (one header line + 9,915 entries):
| counts | regnum |  
|-----:|:-------- |
| 330  | archaea  |
| 7997 | bacteria |
| 1588 | eukaryota |
9915 entries in total

stat_merged.tsv (one header line + 6,521 entries):
| counts | regnum      |  
|-----:  |:----------  |
| 5695   | bacteria    |
| 91     | protists    |
| 59     | plants      |
| 533    | fungi       |
| 49     | metazoa     |
| 94     | vertebrates |
6521 entries in total  

<sup>*</sup>In the annotation from Ensembl bacteria includes also archaea.

---
### main_work
bla, bla, bla
