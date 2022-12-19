# gene_length
## The emergence of eukaryotes signaled by a critical protein-to-gene length ratio

This repository contains the data and programs needed to reproduce the results reported 
in our article: how to obtain the annotations from public repositories, main tables and programs to reproduce the results.  

**The structure of this repository:**  
 - **README.md** guides you all over this repository
 - **main_tables** needed to reproduce the main results  
         -**extra_tables** for the supplementary material and extra information that can be quite helpful (ie. taxonomical ids)
 - **main_work** containing the software needed to reproduce the main results of our work is. This directory also contains:   
         - **main_suppl**, where the programs for the supplementary material are.

---
### Data: the annotations were downloaded from public repositories:

#### Proteins
[Reference proteomes](https://www.uniprot.org/proteomes/?query=*&fil=reference%3Ayes)
were downloaded from [Uniprot](https://www.uniprot.org/). 
Each proteome has a unique Uniprot-identifier (UPID). Also, 
a [description](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/README) 
of the proteomes, as well as a [table](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/README) 
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
our_mnt_dir is the local directory where the data were downloaded.

#### Protein coding genes 
The gene annotations were obtained from different Ensembl's webservers 
for [prokaryotes (archaea, bacteria)](https://bacteria.ensembl.org), [protists](https://protists.ensembl.org), [plants](https://plants.ensembl.org), [fungi](https://fungi.ensembl.org), [metazoa](https://metazoa.ensembl.org), 
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
repositories. For instance, for _Homo sapiens_: 
```
Homo_sapiens.GRCh38.98.gtf.gz @
our_mnt_dir + data/compressed/ + "ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/"
```
our_mnt_dir is, as above, the local directory where the data were downloaded.

##### Taxonomy ids of the different species annotated in Ensembl
The taxonomy id of each species has been downloaded from Ensembl for each division: https://ftp.ensembl.org/pub/release-98/species_EnsemblVertebrates.txt
[prokaryotes (archaea, bacteria)](http://ftp.ensemblgenomes.org/pub/bacteria/release-49/species_EnsemblBacteria.txt),
[protists](http://ftp.ensemblgenomes.org/pub/protists/release-49/species_EnsemblProtists.txt), [plants](http://ftp.ensemblgenomes.org/pub/plants/release-49/species_EnsemblPlants.txt),
[fungi](http://ftp.ensemblgenomes.org/pub/fungi/release-49/species_EnsemblFungi.txt), [metazoa](http://ftp.ensemblgenomes.org/pub/metazoa/release-49/species_EnsemblMetazoa.txt), 
[vertebrates](https://ftp.ensembl.org/pub/release-98/species_EnsemblVertebrates.txt).  

---
### main_tables
For protein coding genes, proteins, and the intersection set between them (merged). The files are provided in standard [tab-separated values](https://en.wikipedia.org/wiki/Tab-separated_values) (*.tsv).
- stat_protCodGenes.tsv (one header line + 33,629 entries)
- stat_proteins.tsv (one header line + 9,915 entries)
- stat_merged (one header line + 6,521 entries)

#### **Number of entries per taxonomical division:**  
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

<sup>*</sup>In the annotation from Ensembl Bacteria includes also Archaea.

stat_proteins.tsv (one header line + 9,915 entries):
| counts | regnum |  
|-----:|:-------- |
| 330  | archaea  |
| 7997 | bacteria |
| 1588 | eukaryota<sup>*</sup> |
9915 entries in total

<sup>*</sup>In the annotations from Uniprot, Eukaryota includes all the clades described above: protists, plants, fungi, metazoa, vertebrates.  

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

### extra_tables
- species_Ensembl.tsv. For each division, the file containing the taxonomy ids of the different species annotated in Ensembl, [see above](./README.md#taxonomy-ids-of-the-different-species-annotated-in-ensembl), has been downloaded. The files for the different divisions have been concatenated into species_Ensembl.tsv, maintaining only the first header. Finally the file has been slimmed-down reducing its columns to cols 1, 2, and 4; that is, colloquial name of the species, species name, and taxonomy_id.

---
### main_work
- protCodGenes_lognormDist.ipynb and proteins_lognormDist.ipynb: the distributions of the lengths of the protein coding genes and proteins, respectively. That is Fig.1 (also S1, S2, and S7)
- protCodGenes_taylorLaw.ipynb and proteins_taylorLaw.ipynb: the observed Taylor law in the lengths distributions for the different species (variance vs mean in $log{_10}$ representation). That is, Fig. 2
