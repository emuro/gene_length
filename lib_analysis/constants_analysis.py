# python3
# ################################################################## #
# constants_analysis.py (C) Jan-2021 Mainz.
# Author: Enrique M. Muro
# ################################################################## #
#
# ------------------------------------------------------------------------
# Project: geneLength
#
# Purpose: Constant to be used in the project
#          This is for the sake of the analysis (It is simplified from constants.py)
# Important:
# Change: DIVISION. There rest adapt to this selection
# ################################################################## #
from datetime import date


DIVISION = "protist"  # "viruses", "bacteria", "protist", "fungi", "plants", "metazoa", "vertebrates"
# ################################################################## #

ENSEMBL_DIVISION = DIVISION
DB = "ensemblgenomes"   # ensembl for vertebrates
                        # ensemblgenomes for bacteria,fungi,protist,plants,metazoa and viruses
ENSEMBL_VERSION = "49"  # these for "bacteria", "protist", "plants", and "metazoa"

if ENSEMBL_DIVISION == "vertebrates":
    ENSEMBL_DIVISION = ""  # division should be here and empty string for the sake of the code
    DB = "ensembl"
    ENSEMBL_VERSION = "98"  # just for vertebrates
elif ENSEMBL_DIVISION == "viruses":
    ENSEMBL_VERSION = "101"  # just for viruses

BOOL_CHECK_SOME_SPECIES = 0  # for debugging
ANNOTATION_NAMES_TO_CHECK = ["Loxodonta_africana",
                             ""]  # ["Homo_sapiens", "Danio_rerio"] ["Homo_sapiens"]

CHECK_GENE_TYPE = ["protein_coding"]   # "[transcribed_processed_pseudogene", "protein_coding", "lncRNA"]


BASE_PATH = "/Volumes/Wes/"  # HOME_PATH; Wes; BIRD
PROJECT = "geneLength"
OUT_DATA_LOCAL_PATH_ROOT = BASE_PATH + "results/" + PROJECT + "/"
OUTPUT_INPUT_FILES_PATH = OUT_DATA_LOCAL_PATH_ROOT + "outputInputFiles/"  # or PYCHARM_PROJECT_PATH

# BIOTYPE_DISTR_PATH_NAME = "biotype_distribution"
BIOTYPE_DISTR_PATH_NAME = "biotype_distribution"
GENES_PATH_NAME = "genes"
DATA_PATH_NAME  = "data/compressed"
SOME_STATISTICS_PATH_NAME = "analysis/some_statistics/stat_description"
SOME_STATISTICS_PROTEINS_PATH_NAME = "analysis/some_statistics/stat_description/proteins/"
SOME_STATISTICS_BIOTYPES_PER_DIVISION_PATH_NAME = "analysis/some_statistics/biotypes_per_division"

# MODEL_ORGANISMS
MODEL_ORGANISMS_PROTEINS_PATH_NAME = "some_tables/model_organisms/proteins/"
MODEL_ORGANISMS_PROTEINS_FILE = "uniprot_well_annotated_organisms.tsv"
REFERENCE_PROTEOMES_PATH_NAME = "some_tables/reference_proteomes/"
REFERENCE_PROTEOMES_FILE = "reference_proteomes_table_28.5.2021.tsv"

# TAXID
ENSEMBL_TAXID_FILE_NAME = "some_tables/species_Ensembl_taxid/species_Ensembl.tsv"

# LOG FILES
############
BIOTYPE_DISTR_OUTPUT_INPUT_PATH = OUTPUT_INPUT_FILES_PATH + BIOTYPE_DISTR_PATH_NAME + "/"
LOG_BIOTYPE_DISTRIBUTION_SPECIES_FILE = BIOTYPE_DISTR_OUTPUT_INPUT_PATH + "biotype_distribution_of_species_" +\
                                   DIVISION + "_" + DB + "_" + ENSEMBL_VERSION + ".tsv"

LOG_GENES_FILE = OUTPUT_INPUT_FILES_PATH + GENES_PATH_NAME + "/" + "division_biotype_species." + DIVISION + "." \
                 + DB + "." + str(ENSEMBL_VERSION) + ".tsv"

LOG_SOME_STATISTICS_FILE = OUTPUT_INPUT_FILES_PATH + SOME_STATISTICS_PATH_NAME + "/" + \
                           "stat_description." + DIVISION + "." + DB + "." + str(ENSEMBL_VERSION) + ".tsv"
LOG_SOME_STATISTICS_TAXID_MERGED_FILE = OUTPUT_INPUT_FILES_PATH + SOME_STATISTICS_PATH_NAME + "/taxid_merged/"\
                           "stat_description." + "taxid_merged.ensembl_and_ref_proteome" + ".tsv"

LOG_SOME_STATISTICS_PROTEIN_FILE = OUTPUT_INPUT_FILES_PATH + SOME_STATISTICS_PROTEINS_PATH_NAME + "/" + \
                           "stat_description." + "protein.uniprot_model_organism" + ".tsv"
LOG_SOME_STATISTICS_REFERENCE_PROTEOME_FILE = OUTPUT_INPUT_FILES_PATH + SOME_STATISTICS_PROTEINS_PATH_NAME + "/" + \
                           "stat_description." + "protein.uniprot_reference_proteome" + ".tsv"


LOG_SOME_STATISTICS_BIOTYPES_PER_DIVISION_FILE = \
    OUTPUT_INPUT_FILES_PATH + SOME_STATISTICS_BIOTYPES_PER_DIVISION_PATH_NAME  + "/" + \
    "biotype_per_division." + str(date.today()) + ".tsv"

