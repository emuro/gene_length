# GeneLength project

## Python3 code (MacBook): 

geneLenght project in the geneLength project dir at enriquem.muro@mac:~/PycharmProjects/geneLength/
general project configuration:
* constants.py: define the constants for the configuration @geneLength/lib/
----

### some_statistics (genes)
@some_statistics:some_statistics.py
Given a division (vertebrates, metazoa, ...), it statistically describes the gene length of each division:biotype:species.  
Note: at the moment is done **only for protein coding**.

It needs:
- The genes of every division:biotype:species (see dir get_genes)


Get the genes for each species (one file per division:biotype:species) 

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
----






