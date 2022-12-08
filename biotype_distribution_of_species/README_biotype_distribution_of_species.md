# GeneLength project

## Python3 code (MacBook): 

geneLenght project in the geneLength project dir at enriquem.muro@mac:~/PycharmProjects/geneLength/
general project configuration:
* constants.py (or constants__arcturus.py): define the constants for the configuration @geneLength/lib/
----

### biotype_distribution_of_species
@ biotype_distribution_of_species:main_biotype_distribution_of_species.py

Given the DIVISION, It gets for each of its species the distribution of the gene types (aka biotypes). 
Annotates the biotype distribution as a line in an specific file for the DIVISION and 
saves the histogram (png plot) for each division:species.

**main_biotype_distribution_of_species.py:** 

* from lib import constants as c  # constants__arcturus
* from lib import files
* from lib import use_pyensembl as pyens
* from lib import plot
* from lib import species as sp

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







