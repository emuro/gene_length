# python3
# ################################################################## #
# main_some_statistics.py (C) March-April-2021 Mainz.
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
import numpy as np
from time import time
import sys
import pathlib

local_biotype = "protein_coding"


def get_biotype_distribution_species(biotype_distribution_file_name):
    """input
    param: biotype_distribution_file_name
    Get the species info and the biotype distribution
    NOTE: the species from these file has been previously indexed (see main_index_ensembl_gtfFiles.py)
    and the indexed file was accessed to obtain the biotype distribution of each species
    (see main_biotype_distribution_of_species andm main_index_ensembl_gtfFiles.py)
    """
    df_biotype_distr = files.get_df_from_tsv_file(biotype_distribution_file_name)
    for i in range(len(df_biotype_distr)):  # Fix typo: iterate through each species and replace "//" by "/"
        df_biotype_distr.loc[i, "trunk_biotype_distribution_path"]= \
            df_biotype_distr.loc[i, "trunk_biotype_distribution_path"].replace("//", "/")
    if 0:
        pd.options.display.max_columns=None
        print(df_biotype_distr.head(5))
        sys.exit("...jut got df_species with the biotype distribution")
    if 0:
        print(df_biotype_distr.columns)

    return df_biotype_distr


def get_division_biotype_species_filter_and_filter(division_biotype_species_file_name, division="", biotype="",
                                                   species="", trunk_annotation_path=""):
    """
    input,
    division_biotype_species_file_name: the input file name with the division_biotype_species
    output,
    df_division_biotype_species: the data frame that stores where the genes of species file is.

    NOTE: the species from these file has been previously indexed (see main_index_ensembl_gtfFiles.py)
    and the indexed file was accessed to obtain the biotype distribution of each species
    (see main_biotype_distribution_of_species and main_index_ensembl_gtfFiles.py).
    The genes for the biotype of the species of the division has been saved in a file
    (see main_get_genes)
    """
    df_div_bio_spec = files.get_df_from_tsv_file(division_biotype_species_file_name)
    if 0:
        pd.options.display.max_columns=None
        print(df_div_bio_spec.head(5))
        sys.exit("...jut got df_species with the biotype distribution")
    if 0:
        print(df_div_bio_spec.columns)

    # filter by division, biotype and species
    if division:
        df_filtered1 = df_div_bio_spec[df_div_bio_spec["division"]==division]
    else:
        df_filtered1 = df_div_bio_spec

    if biotype:
        df_filtered2 = df_filtered1[df_filtered1["biotype"]==biotype]
    else:
        df_filtered2 = df_filtered1

    if species:
        df_filtered3 = df_filtered2[df_filtered2["species"]==species]
    else:
        df_filtered3 = df_filtered2

    if trunk_annotation_path:
        df_filtered4 = df_filtered3[df_filtered3["trunk_annotation_path"]==trunk_annotation_path]
    else:
        df_filtered4 = df_filtered3

    return df_filtered4


def get_genes_length_description(genes_file_with_path,
                                 species="", assembly="", trunk_genes_path="", genes_file=""):
    """
    input,
    df_genes: length data for each species
    output,
    df_description: description of each species
    """
    df_genes = files.get_df_from_tsv_file(genes_file_with_path)
    df_genes["log10ofLength"] = np.log10(df_genes.length)
    df_genes["logofLength"] = np.log(df_genes.length)
    if 0:
        pd.options.display.max_columns=None
        print(df_genes.head(5))
        sys.exit()

    # Comparisons
    #############
    stat_describe_length = df_genes[["length", "log10ofLength", "logofLength"]].describe()
    log10_var = df_genes["log10ofLength"].var()
    log_var = df_genes["logofLength"].var()
    if 0:
        print(stat_describe_length) # count;mean;std;min;25%;50%;75%;max
        print(stat_describe_length["length"]["std"])
        sys.exit()
    list_n = stat_describe_length.length.to_list()
    list_n.insert(3, df_genes["length"].var())  # count, mean, std,...and now, var
    list_log10 = stat_describe_length.log10ofLength.to_list()
    list_log10.pop(0)
    list_log10.insert(2, df_genes["log10ofLength"].var())  # mean, std,...and now, var
    list_log = stat_describe_length.logofLength.to_list()
    list_log.pop(0)
    list_log.insert(2, df_genes["logofLength"].var())  # mean, std,...and now, var
    new_row = [species, assembly, trunk_genes_path, genes_file] + list_n + list_log10 + list_log
    col_names = ["species", "assembly", "trunk_genes_path", "genes_file",
        "count", "mean", "std", "var", "min", "25perc", "50perc", "75perc", "max",
        "log10_mean", "log10_std", "log10_var", "log10_min", "log10_25perc", "log10_50perc", "log10_75perc", "log10_max",
        "log_mean", "log_std", "log_var", "log_min", "log_25perc", "log_50perc", "log_75perc", "log_max"]

    df_genes_description = pd.DataFrame(dict(zip(col_names, new_row)), columns=col_names, index=[0])  # using zip
    return df_genes_description



def main():
    start_time = time()

    #
    # Get the species info and the biotype distribution for each species
    df_biotype_distr = get_biotype_distribution_species(c.LOG_BIOTYPE_DISTRIBUTION_SPECIES_FILE)

    #
    # Where are the genes of the biotype for the species of the distribution?
    # The genes for the biotype of the species of the division has been saved in a file
    # (see main_get_genes)
    df_div_bio_spec = get_division_biotype_species_filter_and_filter(c.LOG_GENES_FILE,
                                                                     "",                # division
                                                                     local_biotype,     # biotype
                                                                     "",                # species
                                                                     "")                # trunk_genes_path
    if 0:
        print(df_div_bio_spec.columns)
    if 0:
        pd.options.display.max_columns=None
        print(df_div_bio_spec.head(5))
        sys.exit()

    df_div_bio_spec.reset_index(inplace=True)
    # df_div_bio_spec.set_index('trunk_annotation_path', inplace=True)
    if 0:
        pd.options.display.max_columns=None
        print(df_div_bio_spec.head(5))
        sys.exit()
    col_names = ["species", "assembly", "trunk_genes_path", "genes_file",
            "count", "mean", "std", "var", "min", "25perc", "50perc", "75perc", "max",
            "log10_mean", "log10_std", "log10_var", "log10_min", "log10_25perc", "log10_50perc", "log10_75perc", "log10_max",
            "log_mean", "log_std", "log_var", "log_min", "log_25perc", "log_50perc", "log_75perc", "log_max"]
    df_genes_description = pd.DataFrame(columns=col_names)
    for i in range(len(df_div_bio_spec)):  # iterate through each species
        if i % 100 == 0:
            print(str(i) + "\t" + df_div_bio_spec.loc[i, "species"])
            
        genes_file_with_path = c.BASE_PATH + df_div_bio_spec.loc[i, "root_genes_path"] +  \
                                    df_div_bio_spec.loc[i, "trunk_genes_path"] + df_div_bio_spec.loc[i, "genes_file"]
        df_new_entry = get_genes_length_description(genes_file_with_path,
                                                            df_div_bio_spec.loc[i,"species"],
                                                            df_div_bio_spec.loc[i, "assembly"],
                                                            df_div_bio_spec.loc[i, "trunk_genes_path"],
                                                            df_div_bio_spec.loc[i, "genes_file"])
        df_genes_description = df_genes_description.append(df_new_entry, ignore_index=True)
    if 0:
        pd.options.display.max_columns=None
        print(df_genes_description.head(3))
        sys.exit()
    # save_file
    pathlib.Path(c.OUTPUT_INPUT_FILES_PATH + c.SOME_STATISTICS_PATH_NAME).mkdir(parents=True,exist_ok=True)
    df_genes_description.to_csv(c.LOG_SOME_STATISTICS_FILE, sep="\t")
    if 1:
        print(str(df_genes_description.shape[0]) + " species has been processed (some_statistics)")
        total_time = time() - start_time
        print("...total time: %.10f seconds." % total_time)


if __name__ == "__main__":
    main()
