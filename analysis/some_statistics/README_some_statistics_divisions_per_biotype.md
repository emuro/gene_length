# GeneLength project

## Python3 code (MacBook): 

geneLenght project in the geneLength project dir at enriquem.muro@mac:~/PycharmProjects/geneLength/
general project configuration:
* constants.py: define the constants for the configuration @geneLength/lib/
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
**Note:** takes ~0.008 seconds
```file_structure
path: results/geneLength/outputInputFiles/analysis/some_statistics/biotypes_per_division/
file: divisions_per_biotype.2021-04-29.tsv
header: total_divisions_per_biotype     total_species_per_biotype
```






