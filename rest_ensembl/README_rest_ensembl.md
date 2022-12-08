# GeneLength project

----
#### rest_ensembl
It access the [rest.ensembl.org](https://rest.ensembl.org "server") server
* main_rest_ensembl.py: get the vertebrate species from the current version and leave them in @geneLength/outputInputFiles/;  
for instance, for the version 102 (jan.2021), rest_ensembl_102_species.tsv
    * constanst.py (or constanst__arcturus.py)
    * rest_ensembl.py: my interface to the rest ensembl api

output table: it runs for the current version (102, 103,...)

    @ /Volumes/Wes/results/geneLength/outputInputFiles/
    enriquem.muro@MacBookPro-EM outputInputFiles % ls   
    rest_ensembl_102_species.tsv 
    rest_ensembl_103_species.tsv
    rest_ensembl_104_species.tsv
----