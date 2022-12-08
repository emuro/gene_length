# python3
# ################################################################## #
# main_some_statistics_biotypes_per_distribution.py (C) March-April-2021 Mainz.
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

# GLOBAL VARIABLES
##################
DIVISIONS = ["viruses", "bacteria", "protist", "fungi", "plants", "metazoa", "vertebrates"]
DIVISIONS.reverse()
DBS = ["ensemblgenomes", "ensemblgenomes", "ensemblgenomes", "ensemblgenomes", "ensemblgenomes", "ensemblgenomes",
       "ensembl"] # for viruses, bateria, ...respectively
DBS.reverse()
ENSEMBL_VERSIONS = ["101", "49", "49", "49", "49", "49", "98"] # for viruses, bateria, ...respectively
ENSEMBL_VERSIONS.reverse()


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

def get_biotype_list_distribution_from_string(dist_string, minimum=0):
    biotypes_dist = dist_string.split(",")
    df_dist = pd.DataFrame([], columns=['biotype','number'])
    for s in biotypes_dist:
        df_dist.loc[len(df_dist.index)] = s.split(":")

    return df_dist
    # number_of_biotypes = len(df_dist.biotype)
    # return df_dist["biotype"].to_list()

def obtain_df_biotypes_of_division_info(division_name, db, ensembl_version, df_species):
    """ obtain analysis of biotypes of the division
        input:
            -division: the name of the division
            -df_species: (all the species of the division from a biotype distribution annotation file)
                columns=['species', 'assembly', 'division', 'db', 'db_version', 'root_biotype_distribution_path',
                'trunk_biotype_distribution_path', 'min_number_of_genes_per_biotype', 'number_of_biotypes',
                'sorted_biotype_distribution']
        out:
            -out string containing the biotype distribution of the division:
                division\tdb\tensembl_version\ttotal_number_of_species_in_division\t
                sorted_biotype_distribution(biotype1:Nsp1,biotype2:Nsp2,...)
    """
    if 0:
        pd.set_option('display.max_columns',None)
        print(df_species.head(5))
        sys.exit()
    #
    #Get a list of all the biotypes of the division
    count_sp_biotypes = {}
    genes_in_biotypes = {} # I do not show this at the end
    total_number_of_genes = 0
    for i in range(len(df_species)):  # iterate through each species
        df_biotypes_sp = get_biotype_list_distribution_from_string(df_species.loc[i, "sorted_biotype_distribution"])
        for j in range(df_biotypes_sp.shape[0]):
            total_number_of_genes += int(df_biotypes_sp.loc[j, "number"])
            if df_biotypes_sp.loc[j, "biotype"] not in count_sp_biotypes.keys():
                genes_in_biotypes[df_biotypes_sp.loc[j, "biotype"]] = int(df_biotypes_sp.loc[j, "number"])
                count_sp_biotypes[df_biotypes_sp.loc[j,"biotype"]] = 1
            else:
                genes_in_biotypes[df_biotypes_sp.loc[j, "biotype"]] += int(df_biotypes_sp.loc[j, "number"])
                count_sp_biotypes[df_biotypes_sp.loc[j,"biotype"]] += 1
    output = division_name + "\t" + db + "\t" + ensembl_version + "\t" + str(df_species.shape[0]) + "\t"
    for k in sorted(count_sp_biotypes, key=count_sp_biotypes.get, reverse=True):  # by value
        output += k + ":" + str(count_sp_biotypes[k]) + ","

    return output.rstrip(output[-1])  # but the last ,



def main():
    start_time = time()
    output_str = "division\tdb\tensembl_version\ttotal_number_of_species_in_division\t" \
                 "sorted_biotype_distribution(biotype1:Nsp1,biotype2:Nsp2,...)"
    for i in range(len(DIVISIONS)): # foreach division
        division = DIVISIONS[i]
        db = DBS[i]
        ensembl_version = ENSEMBL_VERSIONS[i]
        #
        # Get the species and their repective biotype distributions (for each species)
        df_biotype_distribution = \
            get_biotype_distribution_species(c.OUTPUT_INPUT_FILES_PATH + c.BIOTYPE_DISTR_PATH_NAME + "/" +
                                         "biotype_distribution_of_species_" + division + "_" +
                                         db + "_" + str(ensembl_version) + ".tsv")
        if 0:
            print(len(df_biotype_distribution))
            pd.options.display.max_columns=None
            print(df_biotype_distribution.head(5))
            sys.exit("...jut got df_species with the biotype distribution")

        output_str += "\n" + obtain_df_biotypes_of_division_info(division, db, ensembl_version,
                                                                 df_biotype_distribution)

    # save_file
    pathlib.Path(c.OUTPUT_INPUT_FILES_PATH + c.SOME_STATISTICS_BIOTYPES_PER_DIVISION_PATH_NAME).mkdir(parents=True,exist_ok=True)
    f_h =  open(c.LOG_SOME_STATISTICS_BIOTYPES_PER_DIVISION_FILE, "w")
    f_h.write(output_str + "\n")
    f_h.close()
    if 1:
        total_time = time() - start_time
        print("...total time: %.10f seconds." % total_time)


if __name__ == "__main__":
    main()
