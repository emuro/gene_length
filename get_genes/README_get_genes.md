# GeneLength project

## Python3 code (MacBook): 

geneLenght project in the geneLength project dir at enriquem.muro@mac:~/PycharmProjects/geneLength/
general project configuration:
* constants.py: define the constants for the configuration @geneLength/lib/
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
*  
* 

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






