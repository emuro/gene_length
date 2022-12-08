# GeneLength project

----
## Python3 code (MacBook): 
geneLenght project in the geneLength project dir at enriquem.muro@mac:~/PycharmProjects/geneLength/
general project configuration: 
* constants.py: define the constants for the configuration @geneLength/lib/

----
## GeneLength

----
### GeneLength gene annnotations
Data downloaded from ensembl see the next Google doc [README_geneLength_ensemblAnnotations](https://docs.google.com/document/d/1f89JfDs_CB0tQAVRvYqDM5zdB9zDgNhUYeHJhCeN43s/edit#)

----
##### rest_ensembl
See dir rest_ensembl: README_rest_ensembl.md

----
### index_ensembl_gtfFiles
@ **index_ensembl_gtfFiles:main_index_ensembl_gtfFiles.py**

It imports [pyensembl](https://pypi.org/project/pyensembl/ "python module") module to index the *.gtf.gz (compressed) ensembl annotation
files. pyemsembl indexes \*.gtf for the species, given an ensembl Kingdom: vertebrates, bacteria, fungi, protist,
plants, metazoa or viruses (@ constants.py or constants__arcturus.py).
The files are indexed (*.db) in the corresponding directory (where the annotation file *gtf.gz for that species is)

    IMPORTANT NOTE:
        -install pyensembl
            It was easy to install it in Pycharm. But for downloading the data I did the next:
            (see https://pypi.org/project/pyensembl/)
            * pip install pyensembl
            * pyensembl install --release 98 --species hum

Note: a *.gtf file for human is indexed in ~30sec per species. A file needs to be indexed only once.
* main_index_ensembl_gtfFiles.py: 
get the vertebrate species from the current version and leave them in 
@geneLength/outputInputFiles/; for instance, for vertebrates index_ensembl_gtfFiles_vertebrates_species_ensembl_98.tsv
    * constanst.py (or constanst__arcturus.py from arcturus)
    * use_pyensembl.py: my interface to the ensembl api
    * os_utils.py: here I have my own library based on the os module. Here I coded recursiveSearch_leave_subdirs that 
      find all the final dir-leaves given a root file. 
      

      Once that I indexed, these are the sizes I've got (most of the weigth is because of *db files) 
      du -sh *
        166G	bacteria
         20G	fungi
        8.2G	metazoa
         19G	plants
        5.8G	protist
        1.1M	viruses
--

    There are certain bacterias that the parser was not able to index. They were problematic so I avoid indexing them:
    GTF_FILES_TO_EXCLUDE = ["Saccharopolyspora_erythraea_nrrl_2338_gca_000062885.ASM6288v1.49.gtf.gz",
                            "Bacillus_sp_lk2_gca_001043575.ASM104357v1.49.gtf.gz",
                            "Plasticicumulans_acidivorans_gca_003182095.ASM318209v1.49.gtf.gz",
                            "Helicobacter_acinonychis_str_sheeba_gca_000009305.ASM930v1.49.gtf.gz",
                            "Vibrio_europaeus_gca_001695575.ASM169557v1.49.gtf.gz",
                            "Paeniclostridium_sordellii_gca_900000195.UMC2.49.gtf.gz",
                            "Alkalispirochaeta_sphaeroplastigenens_gca_002916695.ASM291669v1.49.gtf.gz",
                            "Rickettsia_prowazekii_str_madrid_e_gca_000195735.ASM19573v1.49.gtf.gz",
                            "Francisella_tularensis_subsp_tularensis_schu_s4_gca_000008985.ASM898v1.49.gtf.gz",
                            "Stenotrophomonas_maltophilia_gca_001676385.ASM167638v1.49.gtf.gz",
                            "Candidatus_phytoplasma_australiense_gca_000069925.ASM6992v1.49.gtf.gz",
                            "Cupriavidus_necator_h16_gca_000009285.ASM928v2.49.gtf.gz",
                            "Chromobacterium_sp_mwu14_2602_gca_002924365.ASM292436v1.49.gtf.gz",
                            "Yersinia_enterocolitica_subsp_enterocolitica_8081_gca_000009345.ASM934v1.49.gtf.gz",
                            "Bordetella_avium_197n_gca_000070465.ASM7046v1.49.gtf.gz",
                            "Bacillus_cereus_gca_001044905.ASM104490v1.49.gtf.gz",
                            "Bacillus_cereus_gca_001044815.ASM104481v1.49.gtf.gz",
                            "Photobacterium_profundum_ss9_gca_000196255.ASM19625v1.49.gtf.gz",
                            "Chromobacterium_sp_lk11_gca_001043705.ASM104370v1.49.gtf.gz",
                            "Azorhizobium_caulinodans_ors_571_gca_000010525.ASM1052v1.49.gtf.gz"]

output tables for each kingdom: vertebrates, bacteria, fungi, protist, plants, metazoa or viruses 

    @ ~/PycharmProjects/geneLength/outputInputFiles/
    enriquem.muro@MacBookPro-EM outputInputFiles % ls
    index_ensembl_gtfFiles_bacteria_species_ensemblgenomes_49.tsv	index_ensembl_gtfFiles_vertebrates_species_ensembl_98.tsv
    index_ensembl_gtfFiles_fungi_species_ensemblgenomes_49.tsv	index_ensembl_gtfFiles_viruses_species_ensemblgenomes_101.tsv
    index_ensembl_gtfFiles_metazoa_species_ensemblgenomes_49.tsv	index_ensembl_gtfFiles_plants_species_ensemblgenomes_49.tsv
    index_ensembl_gtfFiles_protist_species_ensemblgenomes_49.tsv

----
### biotype_distribution_of_species
@ **biotype_distribution_of_species:main_biotype_distribution_of_species.py**

Given the DIVISION, It gets for each of its species the distribution of the gene types (aka biotypes). 
Annotates the biotype distribution as a line in an specific file for the DIVISION and 
saves the histogram (png plot) for each division:species.

#### Input:
A file with the species of the division that has been already indexed
```file_description
For all divisions,  
@ results/geneLength/outputInputFiles/index_ensembl_gtfFiles_vertebrates_species_ensembl_98.tsv
```

A file with the species gene annotation (indexed)
```file_description
For all species of the division,  
@ data/compressed/ftp.ensembl.org/pub/release-98/gtf/gorilla_gorilla/Gorilla_gorilla.gorGor4.98.gtf.db
```

#### Output:
A file with the gene biotype distribution plain text of the division (ie. vertebrates).
```file_description
For instance, for vertebrates  
@ results/geneLength/outputInputFiles/biotype_distribution/biotype_distribution_of_species_vertebrates_ensembl_98.tsv
```

A plot (*png file) with the gene biotype distribution of each species of the division (ie. vertebrates:).
```file_description
For instance, for homo sapiens 
@ results/geneLength/ftp.ensembl.org/pub/release-98/biotype_distribution/homo_sapiens/biotype_distribution_of_species_homo_sapiens_nan_ensembl_98.png
```

----
### get_genes
@ get_genes:main_get_genes.py  
Given a division (vertebrates, metazoa, ...), it saves in files all the genes (for each biotype for each species of the division).
It needs:
- The annotations file (ensembl) already indexed (see dir index_ensembl_gtfFiles)
- The biotype distribution file for each species (see biotype distribution of species)


**main_get_genes.py:** 

Get the genes for each species (one file per division:biotype:species) 

* constanst.py
* files.py
* use_pyensembl
* species.py (for the Class Species)

**Main output file** (division_biotype_species; one line per biotype of the species of the division). Here it saids where the division:biotype:species:genes are saved 

```file_description
path: /Volumes/Wes/results/geneLength/outputInputFiles/genes/protein_coding/ 
file: division_biotype_species.metazoa.ensemblgenomes.49.tsv
header of file: division;species;assembly;db;db_version;root_annotation_path;trunk_annotation_path;annotation_file;root_biotype_distribution_path;trunk_biotype_distribution_path;biotype_distribution_file;root_genes_path;trunk_genes_path;genes_file;min_number_of_genes_per_biotype;number_of_genes_biotype;number_of_genes_total;number_of_species_in_division
```

For instance loxodonta_africana line of division_biotype_species.vertebrates.ensembl.98.tsv points to:
```file_description
path: /Volumes/Wes/results/geneLength/ftp.ensembl.org/pub/release-98/genes/loxodonta_africana/
file: protein_coding.loxodonta_africana.nan.ensembl.98.tsv
header of file: contig;strand;start;end;biotype;gene_id;gene_name;length;diffLength
```


**pseudocode:**
1. get df_species (from the biotype distribution file for the division)

```file_description
df_species
path: results/geneLength/outputInputFiles/
file: index_ensembl_gtfFiles_vertebrates_species_ensembl_98.tsv
```
- Fix a previous typo in column "trunk_biotype_distribution_path" (replace "//" by "/")


2. prepare the main output file (division_biotype_species) (open/hand files: one line per biotype of the species of the division)

```file_description
path: /Volumes/Wes/results/geneLength/outputInputFiles/genes/protein_coding/ 
file: division_biotype_species.metazoa.ensemblgenomes.49.tsv
header of file: division;species;assembly;db;db_version;root_annotation_path;trunk_annotation_path;annotation_file;root_biotype_distribution_path;trunk_biotype_distribution_path;biotype_distribution_file;root_genes_path;trunk_genes_path;genes_file;min_number_of_genes_per_biotype;number_of_genes_biotype;number_of_genes_total;number_of_species_in_division
```

3. For each species (of the division):
* Instance of Class Species
* Get df_all_genes of the species
* Save the counts of genes in each contig (for all genes of the species together)
```file_description
path: /Volumes/Wes/results/geneLength/ftp.ensembl.org/pub/release-98/genes/loxodonta_africana/
file: count_genes_in_contigs.loxodonta_africana.nan.ensembl.98.tsv (header of file: contig;counts)
```
* For each biotype of the species:
1. Save the species:biotype:genes. The genes are sorter by length and the distance to the closest gene in length is calculated too.

```file_description
path: /Volumes/Wes/results/geneLength/ftp.ensembl.org/pub/release-98/genes/loxodonta_africana/
file: protein_coding.loxodonta_africana.nan.ensembl.98.tsv
header of file: contig;strand;start;end;biotype;gene_id;gene_name;length;diffLength
```

2. Calculate/save the count of genes in each contig (species:biotype)

```file_description	
path: /Volumes/Wes/results/geneLength/ftp.ensembl.org/pub/release-98/genes/loxodonta_africana/
file: count_genes_in_contigs.protein_coding_loxodonta_africana.nan.ensembl.98.tsv (header of file: contig;counts)
```

3. Save the main output file (division_biotype_species) (one line per biotype of the species of the division) 

```file_description
path: /Volumes/Wes/results/geneLength/outputInputFiles/genes/ 
file: division_biotype_species.metazoa.ensemblgenomes.49.tsv
header of file: division;species;assembly;db;db_version;root_annotation_path;trunk_annotation_path;annotation_file;root_biotype_distribution_path;trunk_biotype_distribution_path;biotype_distribution_file;root_genes_path;trunk_genes_path;genes_file;min_number_of_genes_per_biotype;number_of_genes_biotype;number_of_genes_total;number_of_species_in_division
```

----
### some_statistics (genes)
@some_statistics:some_statistics.py  
Given a division (vertebrates, metazoa, ...), it statistically describes the gene length of each division:biotype:species.  

> Important notes: 
> 
> log10 and log are the statistical description of the log10 (or natural log) of the lengths. NOT the log10 
> of the description of linear length.
>
>  **At the moment, it is only calculated for protein coding**.


It needs:
- The genes of every division:biotype:species (see dir get_genes)

* constants_analysis.py
* files.py
* use_pyensembl

**pseudocode:**
1. get df_div_bio_spec (from the biotype distribution file for the division) get the division_biotype_species file, each line indicates where the file with genes are located (division:biotype:species) 

```file_description
path: results/geneLength/outputInputFiles/genes/protein_coding/ 
file: division_biotype_species.metazoa.ensemblgenomes.49.tsv
header of file: division;species;assembly;db;db_version;root_annotation_path;trunk_annotation_path;annotation_file;root_biotype_distribution_path;trunk_biotype_distribution_path;biotype_distribution_file;root_genes_path;trunk_genes_path;genes_file;min_number_of_genes_per_biotype;number_of_genes_biotype;number_of_genes_total;number_of_species_in_division
```

Then filter the file and leave only lines by division, biotype, species <-> trunk_genes_path  
2. For each species (of the division):  
    * Get the genes for the specie (biotype and division) and obtain a statistical description. For instance, for protein coding and loxodonta_africana, the line of division_biotype_species.vertebrates.ensemblgenomes.49.tsv points to:
```file_description
path: results/geneLength/ftp.ensembl.org/pub/release-98/genes/loxodonta_africana/
file: protein_coding.genes.loxodonta_africana.nan.ensembl.98.tsv
header of file: contig;strand;start;end;biotype;gene_id;gene_name;length;diffLength
```
3. Save the **main output file** (division_biotype_species) (one line per species of the division:biotype).  
```file_description
path: results/geneLength/outputInputFiles/analysis/some_statistics/stat_description/
file: stat_description.metazoa.ensemblgenomes.49.tsv
header of file: "species", "assembly", "trunk_genes_path", "genes_file",
                "count", "mean", "std", "var", "min", "25perc", "50perc", "75perc", "max",
                "log10_mean", "log10_std", "log10_var", "log10_min", "log10_25perc", "log10_50perc", "log10_75perc", "log10_max",
                "log_mean", "log_std", "log_var", "log_min", "log_25perc", "log_50perc", "log_75perc", "log_max"
```
| counts | file | lost from indexed files |  
|------: | :-----------------: |  ---: |
| 31943 | stat_description.bacteria.ensemblgenomes.49.tsv | 52 |
| 237   | stat_description.protist.ensemblgenomes.49.tsv  |  0 |
| 1014  | stat_description.fungi.ensemblgenomes.49.tsv    |  0 |  
| 96    | stat_description.plants.ensemblgenomes.49.tsv   |  0 |  
| 115   | stat_description.metazoa.ensemblgenomes.49.tsv  |  0 |   
| 227   | stat_description.vertebrates.ensembl.98.tsv     |  0 |  
| 1     | stat_description.viruses.ensemblgenomes.101.tsv |  0 |    
| 33633 | total |                                           52 |
Note every file has 1 line more (header)

----
### some_statistics_biotype_per_division (genes)
@some_statistics:some_statistics_biotypes_per_division.py  
Given a division ("viruses", "bacteria", "protist", "fungi", "plants", "metazoa"or "vertebrates"), it calculates its gene type (aka biotype) distribution  

It needs:
- The biotype distribution for all the species of each division (biotype_distribution_of_species)

* constanst_analysis.py
* files.py

**pseudocode:**
1. For each division (["viruses", "bacteria", "protist", "fungi", "plants", "metazoa", "vertebrates"])
get the biotype distribution from each species (of the division):

```file_description
path: results/geneLength/outputInputFiles/biotype_distribution/
file: biotype_distribution_of_species_vertebrates_ensembl_98.tsv
```

2. Calculate the biotype distribution for all the species of the division

3. Save the **main output file** biotypes of the division (one line per division).
>**file: biotype_per_division.2021-04-29.tsv**
```file_description
path: results/geneLength/outputInputFiles/analysis/some_statistics/biotypes_per_division/
file: biotype_per_division.2021-04-29.tsv
header of file: division;db;ensembl_version;total_number_of_species_in_division;sorted_biotype_distribution(biotype1:Nsp1,biotype2:Nsp2,...)"
```
**Note:** takes ~5 minutes

----
### some_statistics_divisions_per_biotype (genes)
@some_statistics:some_statistics_divisions_per_biotype.py

It calculates the division distribution for each biotype. It also calculates the total number or species per biotype.
The divisions are vertebrate, protist, fungi, etc. The biotypes are protein_coding, tRNA, pseudogenes, etc. The species 
are homo_sapiens, danio_rerio, etc.

It needs (input file):
the main output file from some_statistics:some_statistics_biotypes_per_division (one line per division). See:

```file_description
path: results/geneLength/outputInputFiles/analysis/some_statistics/biotypes_per_division/
file: biotype_per_division.2021-04-29.tsv
header of file: division;db;ensembl_version;total_number_of_species_in_division;sorted_biotype_distribution(biotype1:Nsp1,biotype2:Nsp2,...)"
```
Needs the next lib/modules:
* constants_analysis.py
* files.py

**Main output file** (divisions per biotype, also species_per_biotype only one line).
> file: **divisions_per_biotype.2021-04-29.tsv**  
``
path: results/geneLength/outputInputFiles/analysis/some_statistics/biotypes_per_division/  
header: total_divisions_per_biotype     total_species_per_biotype
``

**Note:** takes ~0.008 seconds
>
----
### plot_mean_var__gLength (gene length)
@ dir plot: plot_mean_var__geneLength: plot_mean_var__gLength.py  
Purpose,  
1. It calculates the plot (mean, var) of the prot. coding gene length (linear, as annotated by Ensembl) using a log10 visualization.  
2. It calculates and shows the regression line.  
3. It calculates also the count vs mean and count vs var (not anymore).  

#### Inputs:
(I) The statistical description of the gene length (all genes by division, takes all divisions)
> All species:  
>@ results/geneLength/outputInputFiles/analysis/some_statistics/stat_description/  
>>file: stat_description.vertebrates.ensembl.98.tsv (or bacteria,...). 

I count directly in the code. Note output numbers are the same than input numbers (no one lost):

| counts| file | lost from indexed files |  
|-----: | :-----------------: |  ---: |
| 31943 | bacteria    | 52 |
| 237   | protist     |  0 |
| 1014  | fungi       |  0 |  
| 96    | plants      |  0 |  
| 115   | metazoa     |  0 |   
| 227   | vertebrates |  0 |  
| 1     | viruses     |  0 |  

or 

(II) ensembl and uniprot merged by taxID
> Species merged by taxID (Ensembl and proteome):  
>@ results/geneLength/outputInputFiles/analysis/some_statistics/stat_description/taxid_merged/  
>> stat_description.taxid_merged.ensembl_and_ref_proteome.tsv  (6521 species in total)

| counts | genes_division |  
|-----:| :---------- |
| 5695 | bacteria    |
| 91   | protist     |
| 533  | fungi       | 
| 59   | plants      |  
| 49   | metazoa     |   
| 94   | vertebrates | 

| counts | prots_superregnum |  
|-----:|:-------- |
| 227  | archaea  |
| 5468 | bacteria |
| 826  | eukaryota |

| counts | merged_division_superregnum |  
|-----: | :--------- |
| 227   | archaea    |
| 5468 | bacteria    |
| 91   | protist     |
| 533  | fungi       | 
| 59   | plants      |  
| 49   | metazoa     |   
| 94   | vertebrates | 

then, (if desired) 
>**filter** by how well the distribution fit a log-norm: kurtosis, skew

#### Output: 
Important!: **The plot is saved manually (one by one)**. 

Analysis of protein length (mean-variance) with log10 visualization. The number of species (points)
represented are exactly those from the input files, unless that an explicit filter (like kurtois 
is applied).

For instance for (I):
>@ results/geneLength/general_results/mean_var/mean_var__all/Ensembl_vs_uniprot/  
>file: **mean_var.ensembl.archaea_bacteria_protists_fungi_plants_metazoa_plants_vertebrates.png**, etc

For (II):
>@ results/geneLength/general_results/mean_var/mean_var__merged_taxid/Ensembl_vs_uniprot/  
>file: **mean_var.ensembl.archaea_bacteria_protists_fungi_plants_metazoa_plants_vertebrates.png**, etc  

----
----
## ProteinLength

----
### ProteinLength annnotations
Data downloaded directly from [Uniprot](https://www.uniprot.org/)  
Reference proteomes where the UPID from the proteome can been obtained [reference proteomes](https://www.uniprot.org/proteomes/?query=*&fil=reference%3Ayes)  
[README Uniprot: reference proteome](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/README)  
Downloaded annotations @  
data/compressed/ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/ (Bacteria/ Eukaryota/)  

For instance:
Eukaryota/UP000005640/UP000005640_9606.fasta.gz (homo sapiens) 
This file has one protein per protein coding gene (ie. 20,614 prots for homo sapiens).  
See the protein selection criteria:  
https://www.uniprot.org/help/canonical_and_isoforms

| counts | regnum |  
|-----:|:-------- |
| 330  | archaea  |
| 8002 | bacteria |
| 1589  | eukaryota |
| 9939  | viruses |
19860 species in total

----
----
### some_statistics_proteins_model_organisms (protein length)
@ analysis: some_statistics: README_some_statistics_proteins_model_organisms.md

It gets the statitistical description of the length of some model organisms that I selected manually (8 species)

#### Inputs:
* Model orgamisms,  
  path: results/geneLength/outputInputFiles/some_tables/model_organisms/proteins/  
  uniprot_well_annotated_organisms.tsv (from P Mier and some that I add like danio rerio)

* Uniprot annotations,  
  path: data/compressed/ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/ (Bacteria/   Eukaryota/)  
  i.e + Eukaryota/UP000005640/UP000005640_9606.fasta.gz
  
#### Output:
The statistical description of the protein length
>@ results/geneLength/outputInputFiles/analysis/some_statistics/stat_description/proteins/  
file: **stat_description.protein.uniprot_model_organism.tsv** (8 entries) 

----
### plot_mean_var__model_organisms_pLength (protein length).
### Note: This result is not that astonishing for a few species instead of thousands of species
@ plot: plot_mean_variance__proteinLength: plot_mean_var__model_organisms__pLength.py  
Purpose,  
1. It calculates the plot (mean, var) of the protein (model organisms) using a log10 visualization.  
2. It calculates and shows the regression line.  
3. It calculates also the count vs mean and count vs var.  

#### Inputs:
The statistical description of the protein length
@ results/geneLength/outputInputFiles/analysis/some_statistics/stat_description/proteins/  
file: stat_description.protein.uniprot_model_organism.tsv (19854 entries)

#### Output: 
Analysis of protein length (mean-variance) with log10 visualization 
> file: **protein_some_model_organisms.png**  
``@ /Volumes/Wes/results/geneLength/general_results/proteinLength/model_organisms/
``

----
### some_statistics_proteins_reference_proteomes (protein length)
@ dir analysis: some_statistics: README_some_statistics_proteins_reference_proteomes.md

It gets the statitistical description of the length of the reference proteomes: count, mean, std, var, ... 

#### Inputs:
* Reference proteomes,  
  path: results/geneLength/outputInputFiles/some_tables/reference_proteomes/
  file: reference_proteomes_table_28.5.2021.tsv (19854 species: archaea 330, bacteria 7997, eukaryota 1588, viruses 9939)

* Uniprot annotations,  
  path: data/compressed/ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/ (Bacteria/   Eukaryota/)
  * i.e + Eukaryota/UP000005640/UP000005640_9606.fasta.gz

#### Output:
The statistical description of the protein length
> file: **stat_description.protein.uniprot_reference_proteome.tsv**  
`` @ results/geneLength/outputInputFiles/analysis/some_statistics/stat_description/proteins/  
(also 19854 entries + header; superregnum: 330 archaea, 7997 bacteria, 1588 eukaryota, 9939 viruses = 9915 + 9939 viruses)
``
----
### plot_mean_var__reference_proteomes_pLength (protein length)
@ dir plot: plot_mean_variance__proteinLength: plot_mean_var__reference_proteomes__pLength.py.  

Purpose,  
1. It calculates the plot (mean, var) of the protein (reference proteomes) using a log10 visualization.  
2. It calculates and shows the regression line.  
3. It calculates also the count vs mean and count vs var.

#### Output: 
Analysis of protein length (mean-variance) with log10 visualization
> @ results/geneLength/general_results/mean_var/mean_var__all/Ensembl_vs_uniprot/  
>file: **mean_var.ref_proteome.archaea_bacteria_eukaryota.png**  
>Also independently for: archaea, bacteria, (archaea and bacteria), eukaryota.

----
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
### main_gene_length_distribution
@ analysis: length_distribution: main_gene_distribution.py

It gets the species (ensembl) from the statistical description of the gene length file.  
For each species,  
  * it calculates the statistics (normal distribution fit): kurtosis, skew, shapiro_wilkinson, 
  D'agostino_pearson, kolmogorov_smirnov.  
  * it generates and saves a *png with the gene Length distribution and the theoretical normal 
  distribution.
  * it also saves the statistics of the normal distribution fit for all the species.

#### Input:
Statistical description from the different divisions of Ensembl (vertebrates, metazoa,...):
Each entry annotates the gene length mean, var, ... of all the genes of an species.
```file_description
For all species,  
@ results/geneLength/outputInputFiles/analysis/some_statistics/stat_description/  
file: stat_description.*.tsv
```
Just for reference: all__stat_description.ensembl.tsv is a join of stat_description_*_ensembl*.tsv  
it has 33633 species + header;  
division: 31943 EnsemblBacteria; 1014 EnsemblFungi; 237 EnsemblProtists; 96 EnsemblPlants;
115 EnsemblMetazoa; 227 EnsemblVertebrates; 1 EnsemblViruses 

This code builds a df with all the previous files (stat_description.*.tsv) and 
annotates the division. The df has 32948 diff species
Note: it has 32948 and not 32950, because there are some species with the same name between divisions.
>it has 32948 diff species names;  
>division: 31261 EnsemblBacteria; 1013 EnsemblFungi; 237 EnsemblProtists; 96 EnsemblPlants;
>113 EnsemblMetazoa; 227 EnsemblVertebrates; 1 EnsemblViruses

#### Outputs:
A *png file with the gene distribution + theoretical normal for each species.
```file_description
For instance, for homo_sapiens  
@ results/geneLength/ftp.ensembl.org/pub/release-98/genes/homo_sapiens/  
file: protein_coding.geneLength_distrib.homo_sapiens.nan.ensembl.98.png
```

A file with the statistics of the normal distribution fit for all the species from the previous output file
```file_description
For instance, for homo_sapiens  
@ results/geneLength/outputInputFiles/analysis/some_statistics/stat_description/  
file: all.logNormStat.description.ensembl.tsv
```
For the original 33633 entries (ensembl species, even if the species name is repeated)

---- 
### main_protein_length_distribution
@ length_distribution: main_protein_distribution.py

It gets the species (Uniprot reference proteomes) from the statistical description of the protein lengths files.    
For each species,  
* it calculates the statistics (normal distribution fit): kurtosis, skew, shapiro_wilkinson, D'agostino_pearson, kolmogorov_smirnov
* it creates and stores  a *png with the protein length distribution and theoretical normal distribution
* it saves the statistics of the normal distribution fit for all the species.

#### Input:
Statistical description from the different reference proteomes (19854):
330 archaea, 7997 bacteria, 1588 eukaryota, 9939 viruses
```file_description
For all species,  
@ results/geneLength/outputInputFiles/analysis/some_statistics/stat_description/proteins/  
file: stat_description.protein.uniprot_reference_proteome.tsv
```

#### Outputs:
A *.png file with the gene distribution + theoretical normal for each species got in the input.
```file_description
For instance, for homo_sapiens  
@ results/geneLength/ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/  
file: UP000005640_9606.length.png
```
Note: I am not plotting viruses because they tend not to have enough proteins  
330 archaea, 7997 bacteria, 1588 eukaryota
The plots are calculated in a little bit more than 3 hours

A file with the statistics of the normal distribution fit for all the species from the previous output file
```file_description
For instance, for homo_sapiens  
@ results/geneLength/outputInputFiles/analysis/some_statistics/stat_description/proteins/  
file: all.logNormStat.description.reference_proteome.tsv (9915 entries + header)
```
Note: I am not calculating for viruses because they tend not to have enough proteins  
330 archaea, 7997 bacteria, 1588 eukaryota (total: 9915) 
This table without plots is calculated in about 6 minutes
---- 