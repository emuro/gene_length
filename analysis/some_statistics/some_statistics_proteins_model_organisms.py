# python3
# ################################################################## #
# main_some_statistics_proteins.py (C) March-April-2021 Mainz.
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
import gzip
from Bio import SeqIO
import numpy as np
from time import time
import sys
import pathlib


def get_prot_length_description(df_prot,
                                 species="", uniprot_fasta_file=""):
    """
    input,
    df_genes: length data for each species
    output,
    df_description: description of each species
    """

    df_prot["log10ofLength"] = np.log10(df_prot.length)
    df_prot["logofLength"] = np.log(df_prot.length)
    if 0:
        pd.options.display.max_columns=None
        print(df_prot.head(5))
        sys.exit()

    # Comparisons
    #############
    stat_describe_length = df_prot[["length", "log10ofLength", "logofLength"]].describe()
    log10_var = df_prot["log10ofLength"].var()
    log_var = df_prot["logofLength"].var()
    if 0:
        print(stat_describe_length) # count;mean;std;min;25%;50%;75%;max
        print(stat_describe_length["length"]["std"])
        sys.exit()
    list_n = stat_describe_length.length.to_list()
    list_n.insert(3, df_prot["length"].var())  # count, mean, std,...and now, var
    list_log10 = stat_describe_length.log10ofLength.to_list()
    list_log10.pop(0)
    list_log10.insert(2, df_prot["log10ofLength"].var())  # mean, std,...and now, var
    list_log = stat_describe_length.logofLength.to_list()
    list_log.pop(0)
    list_log.insert(2, df_prot["logofLength"].var())  # mean, std,...and now, var
    new_row = [species, uniprot_fasta_file] + list_n + list_log10 + list_log
    col_names = ["species", "uniprot_fasta_file",
        "count", "mean", "std", "var", "min", "25perc", "50perc", "75perc", "max",
        "log10_mean", "log10_std", "log10_var", "log10_min", "log10_25perc", "log10_50perc", "log10_75perc", "log10_max",
        "log_mean", "log_std", "log_var", "log_min", "log_25perc", "log_50perc", "log_75perc", "log_max"]

    df_genes_description = pd.DataFrame(dict(zip(col_names, new_row)), columns=col_names, index=[0])  # using zip
    return df_genes_description



def main():
    start_time = time()

    #
    # Get the species info and the biotype distribution for each species
    model_org_file = c.OUTPUT_INPUT_FILES_PATH + c.MODEL_ORGANISMS_PROTEINS_PATH_NAME + c.MODEL_ORGANISMS_PROTEINS_FILE
    l_species = files.get_list_from_column_from_tsv_file(model_org_file, "species")
    if 1:
        print(l_species)
    #
    # model organims data
    df_model_org = pd.read_csv(model_org_file, sep="\t")
    if 0:
        pd.options.display.max_columns=None
        print(df_model_org.head(5))
        sys.exit()

    col_names = ["species", "uniprot_fasta_file",
            "count", "mean", "std", "var", "min", "25perc", "50perc", "75perc", "max",
            "log10_mean", "log10_std", "log10_var", "log10_min", "log10_25perc", "log10_50perc", "log10_75perc", "log10_max",
            "log_mean", "log_std", "log_var", "log_min", "log_25perc", "log_50perc", "log_75perc", "log_max"]
    df_prot_stat_desc = pd.DataFrame(columns=col_names)

    for i in range(len(df_model_org)):  # iterate through each species
        if 1:
            print(str(i) + "\t" + df_model_org.loc[i, "species"])
        uniprot_file_with_path = c.BASE_PATH + c.DATA_PATH_NAME + "/" + df_model_org.loc[i, "uniprot_fasta_file"]
        print(uniprot_file_with_path)

        l_prot_len = []
        with gzip.open(uniprot_file_with_path, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                l_prot_len.append(len(record.seq))
            df_prot_len = pd.DataFrame(l_prot_len , columns=['length'])
            df_new_entry_stat_desc = get_prot_length_description(df_prot_len, df_model_org.loc[i, "species"],
                                                           df_model_org.loc[i, "uniprot_fasta_file"])
            if 0:
                pd.options.display.max_columns=None
                print(df_new_entry_stat_desc.head(3))
                sys.exit()

        df_prot_stat_desc = df_prot_stat_desc.append(df_new_entry_stat_desc, ignore_index=True)
        if 0:
            pd.options.display.max_columns=None
            print(df_prot_stat_desc.head(9))
            sys.exit()
    # save_file
    pathlib.Path(c.OUTPUT_INPUT_FILES_PATH + c.SOME_STATISTICS_PROTEINS_PATH_NAME).mkdir(parents=True,exist_ok=True)
    df_prot_stat_desc.to_csv(c.LOG_SOME_STATISTICS_PROTEIN_FILE, sep="\t", index=False)
    if 1:
        print(str(df_prot_stat_desc.shape[0]) + " species has been processed (some_statistics)")
        total_time = time() - start_time
        print("...total time: %.10f seconds." % total_time)


if __name__ == "__main__":
    main()
