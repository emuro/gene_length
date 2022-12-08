# python3
# ################################################################## #
# constants.py (C) Jan-2021 Mainz.
# Author: Enrique M. Muro
# ################################################################## #
#
# --------------------------------------------------------------------
# Project: geneLength
#
# Purpose: Constant to be used in the project
#
# Important:
# Change: DIVISION. There rest adapt to this selection
#         GENE_TYPE.
# ################################################################## #
from os.path import expanduser  # to get the home path dir

DIVISION = "metazoa"  # "viruses", "bacteria", "protist", "fungi", "plants", "metazoa", "vertebrates"
GENE_TYPE = "protein_coding"   # "transcribed_processed_pseudogene", "protein_coding", "lncRNA"
# ################################################################## #

ENSEMBL_DIVISION = DIVISION
DB = "ensemblgenomes"   # ensembl for vertebrates
                        # ensemblgenomes for bacteria,fungi,protist,plants,metazoa and viruses
ENSEMBL_VERSION = "49"  # these for "bacteria", "protist", "plants", and "metazoa"

if ENSEMBL_DIVISION == "vertebrates":
    ENSEMBL_DIVISION = ""  # division should be here and empty string for the sake of the code
    DB = "ensembl"
    ENSEMBL_VERSION  = "98"  # just for vertebrates
elif ENSEMBL_DIVISION == "viruses":
    ENSEMBL_VERSION = "101"  # just for viruses

BOOL_CHECK_SOME_SPECIES = 0  # for debugging
ANNOTATION_NAMES_TO_CHECK = ["Loxodonta_africana", "Homo_sapiens", "Danio_rerio",
                             ""]  # ["Homo_sapiens", "Danio_rerio"] ["Homo_sapiens"]
                                  # "Caenorhabditis_elegans"
                                  # "Danio_rerio" "Saccharomyces_cerevisiae" "Drosophila_melanogaster"
                                  # "Pan_paniscus" "Mola_mola" "Mus_musculus" "Equus_caballus"
                                  # "Anabas_testudineus" "Triticum_aestivum" species
                                  # ["Staphylococcus_aureus_gca_005774675",
                                  #  "Klebsiella_pneumoniae_gca_900507765",
                                  # ""]
CHECK_GENE_TYPE = ["protein_coding"]   # "[transcribed_processed_pseudogene", "protein_coding", "lncRNA"]
ANNOTATION_NAMES_TO_EXCLUDE = ["Sars_cov_2"]

HOME_PATH = expanduser("~") + "/"  # "/Users/enriquem.muro/"
BASE_PATH = "/Volumes/Wes/"  # HOME_PATH; Wes; BIRD
PROJECT = "geneLength"

DB_PATH_ROOT = "data/compressed/"
DB_LOCAL_PATH_ROOT = BASE_PATH + DB_PATH_ROOT

OUT_DATA_PATH_ROOT = "results/" + PROJECT + "/"  # "results/" + PROJECT + "/"
OUT_DATA_LOCAL_PATH_ROOT = BASE_PATH + OUT_DATA_PATH_ROOT

FTP_PATH = "ftp." + DB + ".org/pub/" + ENSEMBL_DIVISION + "/" # i.e ftp.ensembl.org/pub/metazoa/"

DB_LOCAL_PATH = DB_LOCAL_PATH_ROOT + FTP_PATH
OUT_DATA_LOCAL_PATH = OUT_DATA_LOCAL_PATH_ROOT + FTP_PATH
if ENSEMBL_DIVISION == "viruses":
    pass
else:
    DB_LOCAL_PATH += "release-" + ENSEMBL_VERSION + "/"
    OUT_DATA_LOCAL_PATH += "release-" + ENSEMBL_VERSION + "/"

GTF_LOCAL_PATH = DB_LOCAL_PATH + "gtf/"
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
BIOTYPE_DISTR_FILES_TO_EXCLUDE = []
                                  # ["mycobacterium_dioxanotrophicus_gca_002157835",
                                  # "escherichia_coli_gca_001865295", "sphingobium_indicum_b90a_gca_000264945",
                                  # "brevundimonas_sp_gw460_12_10_14_lb2_gca_001636925",
                                  # "xanthomonas_campestris_pv_musacearum_ncppb_4379_gca_000277895",
                                  # "serratia_sp_3acol1_gca_003668775", "pseudomonas_sp_lph1_gca_002037565",
                                  # "sphingorhabdus_sp_m41_gca_001586275",
                                  # "geobacillus_subterraneus_gca_001618685",
                                  # ""]

PYCHARM_PROJECT_PATH = HOME_PATH + "PycharmProjects/" + PROJECT + "/"
OUTPUT_INPUT_FILES_PATH =  OUT_DATA_LOCAL_PATH_ROOT + "outputInputFiles/"  # or PYCHARM_PROJECT_PATH

#COMPARE_PATH_NAME    = "compare2Benford"
#CHI_SQRT_PATH_NAME   = "chiSqrt"
#CORR_PATH_NAME       = "correlation"
BIOTYPE_DISTR_PATH_NAME = "biotype_distribution"
GENES_PATH_NAME = "genes"
##GENE_TYPE_PER_SPECIE = "geneTypesPerSpecie"
##DATA_PATH_NAME       = "data"

# BIOTYPE_DISTR_LOCAL_PATH = OUT_DATA_LOCAL_PATH + BIOTYPE_DISTR_PATH_NAME + "/"
BIOTYPE_DISTR_OUTPUT_INPUT_PATH = OUTPUT_INPUT_FILES_PATH + BIOTYPE_DISTR_PATH_NAME + "/"
GENES_OUTPUT_INPUT_PATH = OUTPUT_INPUT_FILES_PATH + GENES_PATH_NAME + "/" + GENE_TYPE + "/"

MIN_NUM_OF_GENES_PER_BIOTYPE = 0 # 0 to capture every biotype
LOG_ANNOTATION_INDEXED_SPECIES_FILE = OUTPUT_INPUT_FILES_PATH + "index_ensembl_gtfFiles_" + ENSEMBL_DIVISION + \
                                      "_species_" + DB + "_" + ENSEMBL_VERSION + ".tsv"

LOG_BIOTYPE_DISTRIBUTION_SPECIES_FILE = BIOTYPE_DISTR_OUTPUT_INPUT_PATH + "biotype_distribution_of_species_" +\
                                   DIVISION + "_" + DB + "_" + ENSEMBL_VERSION + ".tsv"

if DIVISION == "vertebrates":
    LOG_ANNOTATION_INDEXED_SPECIES_FILE = OUTPUT_INPUT_FILES_PATH + "index_ensembl_gtfFiles_" + DIVISION + \
                                      "_species_" + DB + "_" + ENSEMBL_VERSION + ".tsv"
