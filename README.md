# The emergence of eukaryotes as an evolutionary algorithmic phase transition

This repository contains the data and programs needed to reproduce the results reported in our article. Here, it is also described how to obtain the annotations from public repositories used in this work.  

**README.md** guides you all over this repository. **The structure of this repository is the next:**   
 - **main_tables** needed to reproduce the main figures.  
        - **suppl_tables** for the supplementary material.  
        - **suppl_tables__extra** contains some extra data that can be helpful (ie. taxonomical ids).

 - **main_work** contains the programs needed to reproduce the main results.   
        - **suppl_work**, where the programs for the supplementary material are.  
        - **suppl_work__extra**, where some extra programs that complement the supplementary material are.

- **gl_lib**  contains libs used by programs of this repository
---
### Data: the annotations were downloaded from public repositories

#### Proteins
The [reference proteomes](https://www.uniprot.org/proteomes/?query=*&fil=reference%3Ayes) were downloaded from the Universal Protein Resource ([Uniprot](https://www.uniprot.org/)). Each proteome has a unique Uniprot-identifier (UPID). A [description](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/README) of the proteomes is also provided. It contains a table with information on every proteome: UPIDs, taxonomy_ids, species names, etc. All the reference proteomes were downloaded from [Uniprot FTP repository](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/) on 28.5.2021. Note that viruses were not downloaded and that Uniprot is updated regularly, every eight weeks.  

Then, for each species a fasta file containing its reference proteome was downloaded, preserving the directory structure of the repository. For instance, for _Homo sapiens_ (UPID: UP000005640 and taxonomy id:9606): 
```
UP000005640_9606.fasta.gz @
our_mnt_dir + /data/compressed/ + "ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/"
```
our_mnt_dir is the local directory where the data were downloaded.

#### Protein coding genes 
The protein coding gene annotations were obtained from different webservers provided by Ensembl: [prokaryotes (Archaea, Bacteria)](https://bacteria.ensembl.org), [protists](https://protists.ensembl.org), [plants](https://plants.ensembl.org), [Fungi](http://fungi.ensembl.org/), [invertebrates](https://metazoa.ensembl.org), [vertebrates](https://ensembl.org/index.html). The categorization in groups of organisms is well established by Ensembl. 


| Ensembl ftp site by Kingdom/division                                          | Release            |  
| :---------------------------------------------------------------------------  | :----------------- |  
| [Archaea, Bacteria](http://ftp.ensemblgenomes.org/pub/bacteria/) | ensemblgenomes 49  |  
| [protists](http://ftp.ensemblgenomes.org/pub/protists/)                       | ensemblgenomes 49  |  
| [plants](http://ftp.ensemblgenomes.org/pub/plants/)                           | ensemblgenomes 49  |  
| [Fungi](http://ftp.ensemblgenomes.org/pub/fungi/)                             | ensemblgenomes 49  |  
| [invertebrates](http://ftp.ensemblgenomes.org/pub/metazoa/)                   | ensemblgenomes 49  |  
| [vertebrates (Vertebrata)](https://ftp.ensembl.org/pub/)                      | ensembl 98         |  

The gzip compressed *.gtf.gz (General Transfer Format) gene annotation files were downloaded for the different species preserving the structure of the directories (FTP Ensembl repositories). For instance, for _Homo sapiens_: 
```
Homo_sapiens.GRCh38.98.gtf.gz @
our_mnt_dir + data/compressed/ + "ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/"
```
our_mnt_dir is, as above, the local directory where all the data were downloaded.

##### Taxonomy ids of the different species annotated in Ensembl
The taxonomy id of each species has been downloaded from the corresponding release from Ensembl for each division: [Archaea, Bacteria](http://ftp.ensemblgenomes.org/pub/bacteria/release-49/species_EnsemblBacteria.txt), [protists](http://ftp.ensemblgenomes.org/pub/protists/release-49/species_EnsemblProtists.txt), [plants](http://ftp.ensemblgenomes.org/pub/plants/release-49/species_EnsemblPlants.txt), [Fungi](http://ftp.ensemblgenomes.org/pub/fungi/release-49/species_EnsemblFungi.txt), [invertebrates](http://ftp.ensemblgenomes.org/pub/metazoa/release-49/species_EnsemblMetazoa.txt), [vertebrates](https://ftp.ensembl.org/pub/release-98/species_EnsemblVertebrates.txt).  

##### Genome quality
The data from [NCBI genome](https://www.ncbi.nlm.nih.gov/genome/) was downloaded (20.6.2022) directly from [NCBI genome reports](https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/).

---
#### The lengths of protein coding genes and proteins
The length of any protein coding gene or protein for all the used species can be accessed from our server:   
[https://genford.uv.es:5001/sharing/P79EcUfhE](https://genford.uv.es:5001/sharing/P79EcUfhE)

---
### main_tables
The files for protein coding genes, proteins, and the intersection set between them (merged) are provided in standard [tab-separated values](https://en.wikipedia.org/wiki/Tab-separated_values) (*.tsv):
- stat_protCodGenes.tsv (header line + 33627 entries).  
- stat_proteins.tsv (header line + 9913 entries).  
- stat_merged (header line + 6519 entries).  

Also, a file for the merged set with the mean gene length vs. rho (fraction of nCDS within the protein coding genes). The entries are ordered by ascending mean gene length:  
- rho_vs_gene.dat (three header lines + 6519 entries).

#### **Number of entries per taxonomical division:**
stat_protCodGenes.tsv (header line + 33627 entries):  

| counts | regnum                       |  
|-----:  |:----------                   |
| 31943  | prokaryotes<sup>*</sup>      |
| 237    | protists                     |
| 96     | plants                       |
| 1014   | Fungi                        |
| 115    | invertebrates                |
| 222    | vertebrates                  |
33627 entries in total  

<sup>*</sup>30714 Bacteria and 1229 Archaea.  

stat_proteins.tsv (header line + 9,913 entries):  

| counts | domain |  
|-----:|:-------- |
| 330  | Archaea  |
| 7997 | Bacteria |
| 1586 | Eukaryota<sup>*</sup> |
9913 entries in total

<sup>*</sup>In the annotations from Uniprot, Eukaryota includes: protists (156), plants (184), Fungi (772), invertebrates (226), and vertebrates (248). The 1586 Eukaryotes were classified using the taxonomic hierarchical classification (downloaded on 19.11.2021) provided by Uniprot and based in the NCBI taxonomy database (see [Lineage](https://www.uniprot.org/help/taxonomy)).

stat_merged.tsv (header line + 6519 entries):  

| counts | regnum      |  
|-----:  |:----------  |
| 5468   | Bacteria    |
| 227    | Archaea     |
| 91     | protists    |
| 59     | plants      |
| 533    | Fungi       |
| 49     | invertebrates |
| 92     | vertebrates |
6519 entries in total  

#### suppl_tables
- stat_protCodGenes_ncbiGenomeAssemblyStatus.tsv. Assembly status for the genomes associated to the Ensembl protein coding genes entries. The file is composed by one header and 33627 entries (rows) with 3 columns: species, ensembl_assembly_accession, assembly_status.  

- gene_length_vs_divergence_time.tsv. Average of $< L_{g} >$ and $< log L_{g} >$ of each group of organisms and their divergence time (Mya) obtained from [Timetree](https://timetree.org).

- protCodGenes_averageLg_perGoOrg.txt.  Groups of organisms  with at least 20 species (to compare with proteomes) and the average $< L_{g} >$ of each group in base pairs.

- proteins_averageLp_perGoOrg.txt. Groups of organisms with at least 20 species (to compare with proteomes) and the average $< L_{p} >$ of each group in amino acids.

#### suppl_tables__extra
- species_Ensembl.tsv. The file contains the taxonomy ids of the different species annotated in Ensembl, [see above](https://github.com/emuro/gene_length/blob/main/README.md#taxonomy-ids-of-the-different-species-annotated-in-ensembl). The files for the different divisions have been concatenated into species_Ensembl.tsv, maintaining only the first header. Finally, the file has been slimmed-down reducing its columns to species, species name and taxonomy_id.  

- 480lognormal.dat. Initial seed for the gene growth model: 5000 genes, lognormally distributed (mean=480).

- Homo_sapiens_CDS_nCDS.xlsx. Data needed to compare the length frequency distribution for coding (CDS) and non-coding (nCDS) genetic sequences, see Extended Data Fig. 8.

---
### main_work
- [protCodGenes_lognormDist.ipynb](https://github.com/emuro/gene_length/blob/main/main_work/protCodGenes_lognormDist.ipynb) and [proteins_lognormDist.ipynb](https://github.com/emuro/gene_length/blob/main/main_work/proteins_lognormDist.ipynb): the distributions of the lengths of the protein coding genes (genes hereafter) and proteins respectively. See Fig.1, also Extended Data Figs. 1 and 7.  

- [protCodGenes_taylorLaw.ipynb](https://github.com/emuro/gene_length/blob/main/main_work/protCodGenes_taylorLaw.ipynb) and [proteins_taylorLaw.ipynb](https://github.com/emuro/gene_length/blob/main/main_work/proteins_taylorLaw.ipynb): the observed Taylor law in the distributions of the lengths of genes and proteins (variance vs mean in $log_{10}$ representation) for the different species. See Fig. 2 and Extended Data Fig. 4.  

- [relation_proteins_protCodGenes_lengths.ipynb](https://github.com/emuro/gene_length/blob/main/main_work/relation_proteins_protCodGenes_lengths.ipynb): threshold in the relationship between the mean gene length and the mean protein length for the different species. See Fig. 3 and Extended Data Fig. 9.  

- [rho_nCDS_within_protCodGenes_lengths.ipynb](https://github.com/emuro/gene_length/blob/main/main_work/rho_nCDS_within_protCodGenes_lengths.ipynb): second-order phase transition in the fraction ($\rho$) of non-coding sequences within protein coding genes with the mean gene length as control parameter. See Fig. 4. 
- allowed_states.f. It calculates the allowed states of Fig. 4.
 
#### suppl_work  
-  [mean_vs_time.ipynb](https://github.com/emuro/gene_length/blob/main/main_work/suppl_work/mean_vs_time.ipynb): for the groups of organisms, it is represented the average of the mean gene lengths against their divergence time from LUCA. Similarly, it is displayed, for those groups of organisms, the average of the mean of the gene lengths' logarithm against the evolutionary divergence time from LUCA. That is,  $\overline{\langle L_{g} \rangle}$ (nt) and $\overline{\langle log L_{g} \rangle}$ (nt) vs. divergence time from LUCA (My). See Extended Data Fig. 3.
 
- [protCodGenes__2nd_order_momentum.ipynb](https://github.com/emuro/gene_length/blob/main/main_work/suppl_work/protCodGenes__2nd_order_momentum.ipynb): the observed generalized Taylor law for the protein coding gene length's distributions for the different genomes:  ($\sigma_{g}^{2} + \langle L_{g} \rangle^{2}$) vs $\langle L_{g} \rangle$ in $log_{10}$ representation, where $\langle L_{g}^{2} \rangle$ is the second order momentum. See Extended Data Fig. 4 that complements the main Fig 2.   

- [proteins__2nd_order_momentum.ipynb](https://github.com/emuro/gene_length/blob/main/main_work/suppl_work/proteins__2nd_order_momentum.ipynb): the same for proteins. That is, the observed generalized Taylor law for the protein length's distributions for the different species:  $(\sigma_{p}^{2} + \langle L_{p}\rangle^{2})$ vs $\langle L_{p} \rangle$ in $log_{10}$ representation; where $\langle L_{p}^{2} \rangle$ is the second order momentum. See Extended Data Fig. 4 that complements the main Fig 2.  

- [meanOfLog_logOfMean.ipynb](https://github.com/emuro/gene_length/blob/main/main_work/suppl_work/meanOfLog_logOfMean.ipynb): for the gene length's distributions of the different genomes. Here we compare the mean of the log of the lengths, $\langle log L \rangle$ , and the log of the mean of lengths, $log \langle L \rangle$, in $log_{10}$ representation. It corresponds to the Extended Data Fig. 5.  

 - [average_mean_lengths__order.ipynb](https://github.com/emuro/gene_length/blob/main/main_work/suppl_work/average_mean_lengths__order.ipynb). Extended Data Fig. 6a: order of the average mean gene lengths for the different groups of organisms. Extended Data Fig 6b: same kind of representation for the average protein lengths.  

 - [meanLg_distribution__perGofOrg.ipynb](https://github.com/emuro/gene_length/blob/main/main_work/suppl_work/meanLg_distribution__perGofOrg.ipynb): distribution of the mean gene lengths (Fungi; 1014 genomes). It corresponds to the Extended Data Fig. 9a. Note that the Extended Data Fig. 9b was calculated using code from the main_work section, see [relation_proteins_protCodGenes_lengths.ipynb](https://github.com/emuro/gene_length/blob/main/main_work/relation_proteins_protCodGenes_lengths.ipynb).  

- reliability_fit.ipynb: calculates the log-likelihood that fits the different distributions. See Extended Data Fig. 2.  

- entropy.f: calculates the entropy of the allowed states. See Extended Data Fig. 10.


#### suppl_work__extra  

- gene_growth_simulator.f: example of simulator of gene growth using a multiplicative stochastic factor.

<!---
- [merged_taylorLaw.ipynb](https://github.com/emuro/gene_length/blob/main/main_work/suppl_work__extra/merged_taylorLaw.ipynb). The observed Taylor law in the merged set for the different species for which we have records in both proteins and protein coding genes (variance vs mean in $log_{10}$ representation). This is an extension of Fig. 2.
---!>




