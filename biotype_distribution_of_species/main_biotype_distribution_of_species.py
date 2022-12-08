# python3
# ################################################################## #
# main_biotype_distribution_of_species.py (C) March-2021 Mainz.
# Author: Enrique M. Muro
# ################################################################## #
#
# ------------------------------------------------------------------------
# Project: geneLength
#
# Purpose:
#
# Important:
# Change: DIVISION @ module constants.py
#
# DIVISION = "bacteria"  # "viruses", "bacteria", "protist", "fungi", "plants", "metazoa", "vertebrates"
# There rest should adapt to this selection.
# Also take into account that the minimum number of genes, for protein_coding) is
# MIN_NUM_OF_GENES_PER_BIOTYPE = 0
# ################################################################## #
import sys
sys.path.append('./lib/')
from lib import constants as c  # constants__arcturus
from lib import files
from lib import use_pyensembl as pyens
from lib import plot
from lib import species as sp

#
import pandas as pd
import pathlib
from time import time
import sys
import os


def get_df_biotype_distribution(my_local_species, out_file_name):
    """inputs:
        -instance of class Species

        output:
        -df_biotype_distribution
            A pandas data frame with the biotype distribution of the gene
    """
    try:
        #
        # get ensembl object of gene annotations from the pyemsembl file that
        # I have previously indexed
        data = pyens.index_db_of_geneAnnotation_from_gtfFile(my_local_species.assembly, my_local_species.name,
                                                             my_local_species.get_annotation_file_with_full_path())
        # get the genes of the species
        genes = pyens.retrieve_geneAnnotationObject_fromIndexedData(data)
        if 0:
            print("There are", len(genes), "genes in total (species: raw data)")
            print(genes[2])  # to see the format
            sys.exit()
        # get df annotation of all species (calculate and include the length)
        # and count the species by biotype (gene type)
        df_all_genes = pyens.geneObjectList_to_df_sorted_by_biotype_and_geneLength(genes)
        if 0:
            pd.set_option('display.max_columns', None)
            print(df_all_genes.head(5))
        #
        local_df_sorted_count_biotype = pyens.count_genes_by_biotype_givenDf(df_all_genes,
                                                                             c.MIN_NUM_OF_GENES_PER_BIOTYPE)
        if 0:
            print(local_df_sorted_count_biotype)
    except Exception as e:
        print(__name__, "\t", get_df_biotype_distribution.__name__, ",", sep="")
        print("\t", e)
        files.add_line_to_existing_file(out_file_name, my_species.name)
        return pd.DataFrame() # an empty dataframe

    del(data)
    del(genes)
    del(df_all_genes)
    return local_df_sorted_count_biotype


#
# MAIN code
start_time = time()

# Get the indexed species info.
# NOTE: the species from these file has been previously indexed
# see (main_index_ensembl_gtfFiles.py)
df_species = files.get_df_from_tsv_file(c.LOG_ANNOTATION_INDEXED_SPECIES_FILE)
if 0:
    pd.options.display.max_columns=None
    print(df_species.head(10))
if 0:
    print(df_species.columns)

# check if some species were previously analyzed
flag_rewrite = 1 # to skip running a previously run big process again
pathlib.Path(c.BIOTYPE_DISTR_OUTPUT_INPUT_PATH).mkdir(parents=True, exist_ok=True)  # in case the dir does not exist
LOG_BIOTYPE_DISTRIBUTION_SPECIES_FILE =c.BIOTYPE_DISTR_OUTPUT_INPUT_PATH+"biotype_distribution_of_species_"+ \
                                       c.DIVISION+ "_"+c.DB+ "_"+c.ENSEMBL_VERSION+ ".tsv"
if 0:
    print(LOG_BIOTYPE_DISTRIBUTION_SPECIES_FILE)
    sys.exit()
df_species_from_prev_process = files.get_df_from_tsv_file_or_create_file(LOG_BIOTYPE_DISTRIBUTION_SPECIES_FILE)
if 0:
    print(df_species_from_prev_process.columns)
if 0:
    print(df_species_from_prev_process.shape)
species_already_ckecked = [] # this is the list that will contain the previously analyzed species
if (df_species_from_prev_process.shape[0] >= 1) :
    species_already_ckecked = df_species_from_prev_process["species"].to_list()
    flag_rewrite = 0
    f = open(LOG_BIOTYPE_DISTRIBUTION_SPECIES_FILE, 'a', buffering=1)  # buffering one line

# prepare file: species to avoid
avoid_this_species_file_with_path = c.BIOTYPE_DISTR_OUTPUT_INPUT_PATH + "biotype_distribution_of_species_" +\
                                   "avoid_this_" + c.DIVISION + "_" + c.DB + "_" + c.ENSEMBL_VERSION + ".txt"
if not os.path.exists(avoid_this_species_file_with_path):
    files.creates_an_empty_file(avoid_this_species_file_with_path)
    f_avoid = open(avoid_this_species_file_with_path, "a")
    f_avoid.write("species\n")
    f_avoid.close()

for i in range(len(df_species)):  # iterate through each species
    if 1:
        print(df_species.loc[i, "species"])
    #
    # Filter some species
    if c.BOOL_CHECK_SOME_SPECIES:  # limiting to some species
        if df_species.loc[i, "species"].capitalize() not in c.ANNOTATION_NAMES_TO_CHECK:
            continue
    if df_species.loc[i, "species"] in c.BIOTYPE_DISTR_FILES_TO_EXCLUDE:
        continue
    # check that it was not previously analyzed
    if len(species_already_ckecked) >= 0:
        if df_species.loc[i, "species"] in species_already_ckecked:
            species_already_ckecked.remove(df_species.loc[i, "species"])
            continue

    if i>35000: # por el tema de la memoria
        break

    print(i, "Species to be analyzed: ", df_species.loc[i, "species"])

    # instance of the species Class
    my_species = sp.Species(df_species.loc[i, "species"], df_species.loc[i, "ensembl_division"],
                            df_species.loc[i, "assembly"], c.DB, df_species.loc[i, "ensembl_version"])
    my_species.init_annotation_os(c.BASE_PATH, c.DB_PATH_ROOT,
                            df_species.loc[i, "db_path"].lstrip('/') + "/",
                            df_species.loc[i, "file"], df_species.loc[i, "db_file"])
    my_species.init_biotype_distribution_os(c.BASE_PATH, c.OUT_DATA_PATH_ROOT, "")
    if 0:
        my_species.show_info()
        sys.exit()

    if flag_rewrite:  # prepare file: output biotype distribution
        # save biotype_distribution_of_the_division
        f = open(LOG_BIOTYPE_DISTRIBUTION_SPECIES_FILE, 'a', buffering=1)  # buffering one line
        f.write(my_species.string_basic_info(1) + "\t" + my_species.string_biotype_distribution_paths(1) + "\t" +
                my_species.string_with_the_biotype_distribution(1) + "\n")  # header
        flag_rewrite = 0

    # get the biotype distribution. If it is not possible -> continue
    df_sorted_count_biotype = get_df_biotype_distribution(my_species, avoid_this_species_file_with_path)
    if df_sorted_count_biotype.shape[0] == 0:
        continue
    if 0:
        print(df_sorted_count_biotype)

    # add to the Class (Species)
    my_species.init_biotype_distribution_from_df(df_sorted_count_biotype,c.MIN_NUM_OF_GENES_PER_BIOTYPE)
    if 0:
        my_species.show_info(1)
        sys.exit()
    #
    # plot the gene type distribution
    BOOL_PLOT_GENETYPE_HIST = 1
    if BOOL_PLOT_GENETYPE_HIST:
        if 0:
            print(my_species.get_biotype_distribution_full_path())
            print(my_species.get_file_name_of_biotype_distribution_histogram())
            sys.exit()
        pathlib.Path(my_species.get_biotype_distribution_full_path()).mkdir(parents=True, exist_ok=True)
        plot_title = my_species.db.capitalize() + " " + my_species.db_version + " (" + my_species.name + ")"
        plot.plot_sorted_histogram_of_biotypes(df_sorted_count_biotype,  plot_title,
                                               1,  # 0|1: display|save file
                                               my_species.get_biotype_distribution_full_path() +
                                               my_species.get_file_name_of_biotype_distribution_histogram())
    f.write(my_species.string_basic_info() + "\t" + my_species.string_biotype_distribution_paths() + "\t" +
            my_species.string_with_the_biotype_distribution() + "\n")

    del(my_species)
    del(df_sorted_count_biotype)
    if 1:
        total_time = time() - start_time
        print("...total time: %.10f seconds." % total_time)
    #
    # filter by gene types with more than MIN_NUM_OF_GENES_PER_BIOTYPE
    # in order to process those genes later
    # if 0:
    #    df_all_genes = df_all_genes[(df_all_genes["biotype"].isin(df_sorted_count_biotype["biotype"].to_list()))]
    #    print("\t...after filtering by MIN_NUM_OF_GENES_PER_BIOTYPE (", c.MIN_NUM_OF_GENES_PER_BIOTYPE,
    #          "): there are ", len(df_all_genes), " genes in total (df_all_genes)", sep="")

f.close()  # close file biotype_distribution_of_the_division
