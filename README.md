# gene_length
## The emergence of eukaryotes signaled by a critical protein-to-gene length ratio

This a repository that contains the code needed to reproduce the results reported 
in our article; we start explaining how to obtain the data necessary to produce those 
results.  

---
### The data was downloaded from public repositories:

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

The original directory structure of the FTP Ensembl repositories was preserved.   
For instance, for Homo sapiens (UPID: UP000005640 and taxonomy id:9606): 
```
Homo_sapiens.GRCh38.98.gtf.gz @
our_mount_dir + data/compressed/ + "ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/"
```
where our_mount_dir is the directory/folder where we download the data.


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
The original directory structure of the FTP repository was preserved.   
For instance, for Homo sapiens (UPID: UP000005640 and taxonomy id:9606): 
```
UP000005640_9606.fasta.gz @
our_mount_dir + /data/compressed/ + "ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/"
```



# Añadir esta tabla a las principales. Hace falta una tabla para genes.
| counts | regnum |  
|-----:|:-------- |
| 330  | archaea  |
| 8002 | bacteria |
| 1589 | eukaryota |
19860 species in total



