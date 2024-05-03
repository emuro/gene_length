# ################################################################## #
# constants.py (C) 2023 Mainz.
# Author: Enrique M. Muro
# ################################################################## #
#
# ------------------------------------------------------------------------
# Project: gene_length
#
# Purpose: constants used in the project
# ################################################################## #

import os
import matplotlib
import numpy as np
import sys

# external USB disk
###################
USB_DISK_NAME = "Nubya" # or Wes
system = list(os.uname())[0]
if system == 'Linux':
    base_path_in = "/media/emuro/" + USB_DISK_NAME + "/" # point to the path where your data is
elif system == 'Darwin':
    base_path_in = "/Volumes/" + USB_DISK_NAME + "/" 

# Paths
#######
GIT_PROJECT_PATH        = os.path.dirname(__file__) + "/../" 
MAIN_TABLES_PATH        = GIT_PROJECT_PATH + "main_tables/"
SUPPL_TABLES_PATH       = MAIN_TABLES_PATH + "suppl_tables/"
EXTRA_TABLES_PATH       = MAIN_TABLES_PATH + "suppl_tables__extra/" 
#
BOOL_WORKING_ON = 0
if BOOL_WORKING_ON:
    WORKING_ON_TABLES_PATH  = GIT_PROJECT_PATH + "working_on_tables/"
#
GENES_PROTS_LENGTH_PATH =  base_path_in + "results/geneLength/"

# File names
############
STAT_G_FILE        = MAIN_TABLES_PATH + "stat_protCodGenes.tsv" 
STAT_P_FILE        = MAIN_TABLES_PATH + "stat_proteins.tsv"
STAT_MERGED_FILE   = MAIN_TABLES_PATH + "stat_merged.tsv"
#
if BOOL_WORKING_ON:
    STAT_MERGED_MEDIAN_MODE_FILE = WORKING_ON_TABLES_PATH + "stat_merged_protCodGenes_median_mode.tsv"
    GROUP_VERT = WORKING_ON_TABLES_PATH +  "stat_protCodGenes_with_ncbiGenomeData_FerGroup.tsv"
#
G_NCBI_GENOME_DATA_FILE       = SUPPL_TABLES_PATH + "stat_protCodGenes_ncbiGenomeAssemblyStatus.tsv"
WRONG_ANNOTATIONS_MERGED_FILE = EXTRA_TABLES_PATH + "noisy_stat_merged.tsv"


# colors for plots
################## colors start
COLOR_FOR_DIST = {
    "genes":    matplotlib.colors.to_hex("#76bdfb", keep_alpha=False),
    "proteins": matplotlib.colors.to_hex("#ffab98", keep_alpha=False)
}
#
ORG_KINGDOMS            = ['archaea',  'bacteria', 'eukaryota']
COLOR_KINGDOMS          = ['#D83B01',  '#002050',   '#0078D7']
#
ORG_GROUPS                = ["bacteria", "archaea", "protists", "fungi", "plants",
                             "invertebrates", "vertebrates"]
PyNB_COLOR_ORG_GROUPS     = ['#D83B01',  '#002050', '#FFA500', '#A80000', '#107C10',
                             '#EF008C', '#0078D7']
COLOR_OF = dict([(key, val) for i, (key, val) in enumerate(zip(ORG_GROUPS, PyNB_COLOR_ORG_GROUPS))])
#
ORG_TIME_GROUPS           = ['Archaea', 'Bacteria', 'protists', 'Fungi', 'Viridiplantae', \
        'Arthropoda', 'Actinopterygii', 'Aves', 'Mammalia', 'Primates']
COLOR_ORG_TIME_GROUPS     = [ '#002050', '#D83B01', '#FFA500', '#A80000', '#107C10', '#EF008C', 'powderblue', 'lightskyblue', 'deepskyblue', 'blue']
# 
ORG_MEAN_GROUPS           = ['Archaea', 'Bacteria', \
    'Microsporidia', 'Saccharomycotina', 'Mucuromycota', 'Ascomycota except Sacch.', 'Basidiomycota', \
        'protists', 'Viridiplantae', 'Arthropoda', 'Actinopterygii', 'Aves', 'Mammalia but Primates', 'Primates']
COLOR_ORG_MEAN_GROUPS     = ['#002050', 'slategray', \
    'mistyrose', 'salmon', 'darksalmon', 'coral', 'orangered', \
    '#FFA500', '#107C10', '#EF008C', 'powderblue', 'lightskyblue', 'deepskyblue', 'blue']


# alpha parameter for plots 
# not in use at the moment
################################# alpha start    
ALPHA_OF ={
    "bacteria": 0.1, "archaea": 1.0,
    "protists": 1.0, "fungi": 0.1,
    "plants": 1.0,
    "invertebrates": 1.0, "vertebrates": 1.0
}
ALPHA_ORG_GROUPS = np.array(
    [ALPHA_OF["bacteria"], ALPHA_OF["archaea"],
     ALPHA_OF["protists"], ALPHA_OF["fungi"],
     ALPHA_OF["plants"],
     ALPHA_OF["invertebrates"], ALPHA_OF["vertebrates"]]
)
################################# alpha end




