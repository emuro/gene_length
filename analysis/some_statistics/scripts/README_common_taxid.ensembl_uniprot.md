# GeneLength project

## Python3 code (MacBook): 

geneLenght project in the geneLength project dir at enriquem.muro@mac:~/PycharmProjects/geneLength/

----
### common_taxid.ensembl_uniprot (genes, proteins)
@ analysis: some_statistics: scripts: common_taxid.ensembl_uniprot.py  

>Reduce the species to those whose gene (Ensembl) and protein (Uniprot) annotations  
>can be associated through the taxid for the species name.
>The code provides for the selected species the statistical description of the gene 
>and protein lengths.
> 
>Finally the table is filtered again analyzing the annotations and increasing
>the quality of the results and reducing the number of species. For instance,  
>if the average protein is longer than the gene, etc.

Last filter:
0. From a total of 7671 species: archaea (283), bacteria (6459), fungi (566), metazoa (63), plants (72), protist (114), vertebrates (114).  
1. Filter those species when 3*prots_mean > genes_mean (total: 6707)  
2. Filter when count (genes o prots) < 500 (total: 6685)  
3. Filter when 0.095 > (prots_count/genes_count) or (prots_count/genes_count) > 1.05  (total 6521).   
4. That is, we end up with 6521 species: archaea (227), bacteria (5468), fungi (533), protist (91), plants (59), metazoa (49), vertebrates (94)
Running time: Takes little time (seconds)

Reference_proteomes from Uniprot. See:
* https://www.uniprot.org/proteomes/?query=*&fil=reference%3Ayes  
* ftp://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/reference_proteomes  

#### Inputs:
* gene length Ensembl (statistical description of the protein length: mean, sigma, var, ...)
```file_description
path: results/geneLength/outputInputFiles/analysis/some_statistics/stat_description/
file: all__stat_description.ensembl.tsv (this is a join of stat_description_*_ensembl*.tsv)
```
33633 species + header; 
division: 31943 EnsemblBacteria; 1014 EnsemblFungi; 237 EnsemblProtists; 96 EnsemblPlants; 115 EnsemblMetazoa; 227 EnsemblVertebrates

* species_Ensembl (with tax_id; obtained from Ensembl),
```file_description
path: results/geneLength/outputInputFiles/some_tables/species_Ensembl_taxid/
file: species_Ensembl.tsv (this is a join of species_Ensembl*.txt, obtained by ftp from the corresponding Ensembl release)
```
33021 species + header  
division: 31332 EnsemblBacteria; 1014 EnsemblFungi; 237 EnsemblProtists; 96 EnsemblPlants; 115 EnsemblMetazoa; 227 EnsemblVertebrates

* Reference proteomes (statistical description: protein length: mean, sigma, var, ...)
```file_description
path: results/geneLength/outputInputFiles/analysis/some_statistics/stat_description/proteins/
file: stat_description.protein.uniprot_reference_proteome.tsv 
```
19854 species + header  
superregnum: 330 Archaea; 7,997 Bacteria; 1,588 Eukaryota; 9,939 Viruses =  
= 9915 + 9939 Viruses

#### Output:
The statistical description of the gene length (ensembl) and its associated protein length 
(uniprot; reference proteomes): each species has a different taxid and a different species name.  
```file_description  
@ results/geneLength/outputInputFiles/analysis/taxid_merged/  
file: stat_description.taxid_merged.ensembl_and_ref_proteome.tsv
```
6521 species (unique species name, unique tax id):  
* by gene_division: bacteria (5695), fungi (533), protist (91), plants (59), metazoa (49), vertebrates (94)  
* by prots_superregnum: archaea (227), bacteria (5468), eukaryota (826)  
* by merged_division_superregnum: archaea (227), bacteria (5468), fungi (533), protist (91), plants (59), metazoa (49), vertebrates (94)
   






