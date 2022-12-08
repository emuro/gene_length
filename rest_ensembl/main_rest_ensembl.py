# python3
# ################################################################## #
# main_rest_ensembl.py (C) Jan-2021 Mainz.
# Author: Enrique M. Muro
# ################################################################## #
#
# ------------------------------------------------------------------------
# Project: geneLength
#
# Purpose: get the vertebrate species from the current version and leave them in
# @geneLength/lib/; for instance, for the version 102 (jan.2021)
# it was rest_ensembl_102_species.tsv
#
# Use the next modules:
#   constanst.py (or constanst__arcturus.py)
#   rest_ensembl.py: my interface to the ensembl api
#
# Output file:
#   rest_ensembl_??_species.tsv @ c.OUTPUT_INPUT_FILES_PATH
#   for instance,
#   rest_ensembl_102_species.tsv @ /media/emuro/Wes/results/geneLength/outputInputFiles/ (linux)
#                                @ /Volumes/Wes/results/geneLength/outputInputFiles/ (macbook)
#
# (this is for vertebrates. Also see "https://rest.ensembl.org")
# Future:
# It will be interesting to see if this can be done with bacterias, fungi, etc.
#

import sys
sys.path.append('./lib/')
from lib import rest_ensembl, constants as c  # constants__arcturus


# REST (diff. species)
######################
BOOL_REST_ENSEMBL_SPECIES = 1
if BOOL_REST_ENSEMBL_SPECIES:
    ensembl_release, number_of_species, df_species = rest_ensembl.get_REST_ensembl_species()
    print("\nREST ENSEMBL (rest_ensembl.get_REST_ensembl_species):")
    print("\tThe data corresponds to the release:", ensembl_release, "from ensembl")
    print("\tThere are", number_of_species, "different species in this release")
    print("\tdf_species sorted by taxonId:")
    rest_ensembl_species_file = c.OUTPUT_INPUT_FILES_PATH + \
                                "rest_ensembl_" + str(ensembl_release) + "_species.tsv"
    df_species.to_csv(rest_ensembl_species_file, index = False, sep="\t")
    print("\t", df_species)
    exit()
