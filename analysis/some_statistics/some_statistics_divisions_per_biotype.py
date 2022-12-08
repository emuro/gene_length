# python3
# ################################################################## #
# main_some_statistics_divisions_per_biotype.py (C) March-April-2021 Mainz.
# Author: Enrique M. Muro
# ################################################################## #
#
# ------------------------------------------------------------------------
# Project: geneLength
#
# Purpose:
#
# Important:
# Change:
# ################################################################## #
from lib_analysis import constants_analysis as c
from lib import files
#
import pandas as pd
from time import time
import sys
import pathlib

DATE = "2021-04-29"
local_LOG_SOME_STATISTICS_BIOTYPES_PER_DIVISION_FILE = \
        c.OUTPUT_INPUT_FILES_PATH + c.SOME_STATISTICS_BIOTYPES_PER_DIVISION_PATH_NAME  + "/" + \
        "biotype_per_division." + DATE + ".tsv"
local_LOG_SOME_STATISTICS_DIVISIONS_PER_BIOTYPE_FILE = \
        c.OUTPUT_INPUT_FILES_PATH + c.SOME_STATISTICS_BIOTYPES_PER_DIVISION_PATH_NAME  + "/" + \
        "divisions_per_biotype." + DATE + ".tsv"

def main():
    start_time = time()
    #
    # Get the species and their repective biotype per divisions
    local_LOG_SOME_STATISTICS_BIOTYPES_PER_DIVISION_FILE = \
        c.OUTPUT_INPUT_FILES_PATH + c.SOME_STATISTICS_BIOTYPES_PER_DIVISION_PATH_NAME  + "/" + \
        "biotype_per_division." + DATE + ".tsv"
    df_biotype_per_divisions = files.get_df_from_tsv_file(local_LOG_SOME_STATISTICS_BIOTYPES_PER_DIVISION_FILE)
    if 0:
        print(len(df_biotype_per_divisions))
        pd.options.display.max_columns=None
        print(df_biotype_per_divisions.head(5))
        sys.exit("...jut got df_species with the biotype_per_divisions")

    dict_distr_species = {}
    dict_distr_divisions = {}
    output_species = ""
    output_divisions = ""
    for i in range(df_biotype_per_divisions.shape[0]):
        division = df_biotype_per_divisions.loc[i, "division"]
        distribution = df_biotype_per_divisions.loc[i, "sorted_biotype_distribution(biotype1:Nsp1,biotype2:Nsp2,...)"]
        for s in distribution.split(","):
            [b, n] = s.split(":")
            if b in dict_distr_species: # it was a key? or it is a new one?
                dict_distr_species[b] += int(n)
                dict_distr_divisions[b] += 1
            else:
                dict_distr_species.update({b:int(n)})
                dict_distr_divisions.update({b: 1})
    for k in sorted(dict_distr_divisions, key=dict_distr_divisions.get, reverse=True): # by value
        output_divisions += k + ":" + str(dict_distr_divisions[k]) + ","
        output_species += k + ":" + str(dict_distr_species[k]) + ","
    output = "total_divisions_per_biotype" + "\t" + "total_species_per_biotype" + "\n" + \
             output_divisions.rstrip(output_divisions[-1]) + "\t" + output_species.rstrip(output_species[-1]) + "\n"
    #
    # save_file
    print(output)
    pathlib.Path(c.OUTPUT_INPUT_FILES_PATH + c.SOME_STATISTICS_BIOTYPES_PER_DIVISION_PATH_NAME).mkdir(parents=True,exist_ok=True)
    f_h = open(local_LOG_SOME_STATISTICS_DIVISIONS_PER_BIOTYPE_FILE , "w")
    f_h.write(output)
    f_h.close()
    if 1:
        total_time = time() - start_time
        print("...total time: %.10f seconds." % total_time)


if __name__ == "__main__":
    main()
