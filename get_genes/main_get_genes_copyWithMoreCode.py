# python3
# ################################################################## #
# main_get_genes.py (C) March-April-2021 Mainz.
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
from lib import constants as c
from lib import files
from lib import use_pyensembl as pyens
from lib import species as sp
#
import pandas as pd
import numpy as np
import glob
import pathlib
from time import time
import sys


def get_biotype_list_distribution_from_string(dist_string, minimum=0):
    biotypes_dist = dist_string.split(",")
    df_dist = pd.DataFrame([], columns=['biotype','number'])
    for s in biotypes_dist:
        df_dist.loc[len(df_dist.index)] = s.split(":")
    # number_of_biotypes = len(df_dist.biotype)
    # del df_dist
    return df_dist["biotype"].to_list()

def get_path_species_of_biotype_for_division(base_path="unknown_base_path", root_path="unknown_root_path",
                                             trunk_path="unknown_trunk_path", species="unknown_species"):
    """ i.e:  species_of_biotype.protein_coding.nan.ensembl.98.tsv (nan is for vertebrates)
    """
    f_path = base_path + root_path + trunk_path
    f_path = f_path.replace(species + "/", "")
    f_path = f_path.replace(c.BIOTYPE_DISTR_PATH_NAME, c.GENES_PATH_NAME)
    return f_path

def get_filename_of_species_of_biotype_for_division(biotype="unknown_biotype", division="unknown_division",
                                                    db="unknown_db", db_version="unknown_version"):
    """ i.e:  species_of_biotype.protein_coding.nan.ensembl.98.tsv (nan is for vertebrates)
    """
    return "species_of_biotype." + biotype + "." + str(division) + "." + db + "." + str(db_version) + ".tsv"

def get_filename_of_division_biotype_species(division="unknown_division", db="unknown_db", db_version="unknown_version"):
    """ i.e:  division_biotype_species.vertebrates.ensembl.98.tsv (nan is for vertebrates)
    """
    if str(division) == "nan" and c.DIVISION == "vertebrates":
        return "division_biotype_species." + c.DIVISION + "." + db + "." + str(db_version) + ".tsv"
    else:
        return "division_biotype_species." + str(division) + "." + db + "." + str(db_version) + ".tsv"

def get_header_for_biotypes_of_division_file():
    """ get the header for the biotypes of division files
     """
    return "species\tassembly\tdivision\tdb\tdb_version\troot_genes_path\ttrunk_genes_path\tgenes_file\t" + \
           "min_number_of_genes_per_biotype\tnumber_of_genes\tbiotype\n"

def get_header_for_division_biotype_species_out_file():
    """ get the header for the main output file: division_biotype_species file
    """
    return "division\tbiotype\tspecies\tassembly\tdb\tdb_version\troot_annotation_path\ttrunk_annotation_path\t" \
           "annotation_file\troot_biotype_distribution_path\ttrunk_biotype_distribution_path\t" \
           "biotype_distribution_file\troot_genes_path\ttrunk_genes_path\tgenes_file\tmin_number_of_genes_per_biotype\t" \
           "number_of_genes_biotype\tnumber_of_genes_total\tnumber_of_species_in_division\n"

def obtain_df_biotypes_of_division_info(df_species):
    """ obtain all the biotypes of the division (open/hand files: one per biotype)
        input:
            -df_species: (all the species of the division from a biotype distribution annotation file)
        out:
            -df_biotypes_of_division
                columns=["biotype", "biotype_of_div_file_path", "species_with_biotype_files",
                        "species_with_biotype_file_handler"])
    """
    #Get a list of all the biotypes of the division
    diff_biotypes = []
    for i in range(len(df_species)):  # iterate through each species
        if c.BOOL_CHECK_SOME_SPECIES:  # Filter some species
            if df_species.loc[i, "species"].capitalize() not in c.ANNOTATION_NAMES_TO_CHECK:
                continue
        for i in get_biotype_list_distribution_from_string(df_species.loc[i, "sorted_biotype_distribution"]):
            if i not in diff_biotypes:
                diff_biotypes.append(i)
    diff_biotypes.sort()  # lexicographically

    df_biotype_inf = pd.DataFrame([], columns=["biotypes", "biotype_of_div_file_path",
                                               "biotype_of_div_file", "biotype_of_div_file_handler"])
    for b in diff_biotypes:
        # get path and files
        f_path = get_path_species_of_biotype_for_division(c.BASE_PATH,
                                                          df_species.loc[0, "root_biotype_distribution_path"],
                                                          df_species.loc[0, "trunk_biotype_distribution_path"],
                                                          df_species.loc[0, "species"])
        f_path = f_path.replace("/results/","/testingResults/")  # BORRAR !!!!! solo para testear
        pathlib.Path(f_path).mkdir(parents=True,exist_ok=True)  # just in case it was not previously created
        f = get_filename_of_species_of_biotype_for_division(b, df_species.loc[0, "division"], df_species.loc[0,"db"],
                                                            df_species.loc[0, "db_version"]) # 0 is the first species
        f_h = open(f_path + f, "w")  # open file and add the header
        f_h.write(get_header_for_biotypes_of_division_file())
        new_row = {"biotypes":b, "biotype_of_div_file_path":f_path, "biotype_of_div_file":f,
                   "biotype_of_div_file_handler":f_h}
        df_biotype_inf = df_biotype_inf.append(new_row, ignore_index=True)
        if 0:
            pd.set_option('display.max_columns', None)
            print(df_biotype_inf.head(5))
            sys.exit()
    return df_biotype_inf

def calculate_length_diff_to_closest_gene(df_genes):
    """ Calculate the difference of length to the closest (in length) gene
        input:
            -df_genes: (of the same biotype)
                all the genes are of the same biotype
                already sorted by length
        out:
            -add the next columns to df_genes
                diffLength -> diff of length to the closest (in length) gene
                # log10ofLength -> log (base10) of the length
                # ogofLength -> log (natural) of the length
     """
    l_list = np.array(df_genes['length'].to_list())
    lprev = np.roll(l_list, 1)
    lnext = np.roll(l_list, -1)
    d = lnext - l_list  # distance with the next
    d[-1] = 999999999
    dprev = l_list - lprev  # distance with the prev
    dprev[0] = 999999999
    d[np.less(dprev, d)] = dprev[np.less(dprev, d)]
    df_genes["diffLength"] = d
    # df_genes_of_biotype["log10ofLength"] = np.log10(df_genes_of_biotype.length)
    # df_genes_of_biotype["logofLength"] = np.log(df_genes_of_biotype.length)
    if 0:
        pd.options.display.max_columns = None
        print(df_genes.head(5))


#
# MAIN code
start_time = time()
#
# Get the species info.
# NOTE: the species from these file has been previously indexed (see main_index_ensembl_gtfFiles.py)
# and the indexed file was accessed to obtain the biotype distribution of each species
# (see main_biotype_distribution_of_species andm main_index_ensembl_gtfFiles.py)
df_species = files.get_df_from_tsv_file(c.LOG_BIOTYPE_DISTRIBUTION_SPECIES_FILE)
for i in range(len(df_species)):  # Fix typo: iterate through each species and replace "//" by "/"
    df_species.loc[i,"trunk_biotype_distribution_path"] = \
        df_species.loc[i, "trunk_biotype_distribution_path"].replace("//", "/")
if 0:
    pd.options.display.max_columns = None
    print(df_species.head(5))
    sys.exit("...jut got df_species with the biotype distribution")
if 0:
    print(df_species.columns)

#
# Get all the biotypes of the division (open/hand the corresponding files: one per biotype)
# df columns=["biotypes", "biotype_of_div_file_path", "biotype_of_div_file", "biotype_of_div_file_handler"])
df_biotypes_of_div = obtain_df_biotypes_of_division_info(df_species)
if 0:
    pd.set_option('display.max_columns', None)
    print(df_biotypes_of_div.head(5))
    sys.exit()

#
# prepare the main output file (division_biotype_species)
division_biotype_species_out_file = get_filename_of_division_biotype_species(df_species.loc[i, "division"],
                                                                             df_species.loc[i, "db"],
                                                                             df_species.loc[i, "db_version"])
pathlib.Path(c.GENES_OUTPUT_INPUT_PATH ).mkdir(parents=True,exist_ok=True)  # just in case it was not previously created
division_biotype_species_f_h = open(c.GENES_OUTPUT_INPUT_PATH + division_biotype_species_out_file, "w")  # open file
division_biotype_species_f_h.write(get_header_for_division_biotype_species_out_file())
if 1:
    print(c.GENES_OUTPUT_INPUT_PATH + "\t" + division_biotype_species_out_file)


for i in range(len(df_species)):  # iterate through each species
    if 1:
        print(df_species.loc[i, "species"])
    #
    # Filter some species
    if c.BOOL_CHECK_SOME_SPECIES:  # limiting to some species
        if df_species.loc[i, "species"].capitalize() not in c.ANNOTATION_NAMES_TO_CHECK:
            continue

    # instance of the Class species
    my_species = sp.Species(df_species.loc[i, "species"], df_species.loc[i, "division"],
                            df_species.loc[i, "assembly"], c.DB, df_species.loc[i, "db_version"])
    db_path = df_species.loc[i, "trunk_biotype_distribution_path"].replace(c.BIOTYPE_DISTR_PATH_NAME, "gtf", 1)
    full_db_path = c.BASE_PATH + c.DB_PATH_ROOT + db_path
    db_file = glob.glob(full_db_path + "*.db")[0].replace(full_db_path, "")

    my_species.init_annotation_os(c.BASE_PATH, c.DB_PATH_ROOT,
                                  db_path, db_file.replace(".db", ".gz"), db_file)
    my_species.init_biotype_distribution_os(c.BASE_PATH, c.OUT_DATA_PATH_ROOT,
                                            df_species.loc[i, "trunk_biotype_distribution_path"])
    my_species.init_biotype_distribution_from_string(df_species.loc[i, "sorted_biotype_distribution"],
                                                     df_species.loc[i, "min_number_of_genes_per_biotype"])
    my_species.init_genes(c.BASE_PATH, c.OUT_DATA_PATH_ROOT,
        df_species.loc[i, "trunk_biotype_distribution_path"].replace(c.BIOTYPE_DISTR_PATH_NAME, c.GENES_PATH_NAME, 1))
    if 0:
        my_species.show_info(1)
        sys.exit()

    # get ensembl object of gene annotations from the pyemsembl file that
    # I have previously indexed
    data = pyens.index_db_of_geneAnnotation_from_gtfFile(my_species.assembly, my_species.name,
                                                         my_species.get_annotation_file_with_full_path())
    # get the genes of the species
    genes = pyens.retrieve_geneAnnotationObject_fromIndexedData(data)
    del data
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
        sys.exit()
    #
    count_genes_in_contig_file = my_species.get_genes_of_species__full_path() + \
                                 my_species.get_species_file_name_of_genes_count_foreach_contig()
    s_all = df_all_genes["contig"].value_counts().reset_index()  # sort by frequencies by default
    s_all.columns = ['contig', 'counts']
    pathlib.Path(my_species.get_genes_of_species__full_path()).mkdir(parents=True, exist_ok=True)
    s_all.to_csv(count_genes_in_contig_file, sep="\t", index=True)

    # For each biotype
    gp = df_all_genes.groupby('biotype')

    for g in gp.groups:
        # save the genes
        genes_of_biotype_file = my_species.get_genes_of_species__full_path() + \
                                my_species.get_file_name_of_genes_of_species_and_biotype(str(g))
        df_genes_of_biotype = pd.DataFrame(gp.get_group(g))
        df_genes_of_biotype.reset_index(drop=True, inplace=True)
        number_of_genes = df_genes_of_biotype.shape[0]
        calculate_length_diff_to_closest_gene(df_genes_of_biotype)
        df_genes_of_biotype.to_csv(genes_of_biotype_file, sep="\t", index=True)

        # save the count of genes in each contig
        s = df_genes_of_biotype["contig"].value_counts().reset_index()  # sort by frequencies by default
        s.columns = ['contig', 'counts']
        biotype_genes_in_contig_file = my_species.get_genes_of_species__full_path() + \
                                       my_species.get_biotype_file_name_of_genes_count_foreach_contig(str(g))
        s.to_csv(biotype_genes_in_contig_file, sep="\t", index=True)

        # save entry of species in division::biotype file
        row_df = df_biotypes_of_div.loc[df_biotypes_of_div["biotypes"]==str(g),].iloc[0] # iloc: there is only 1 row
        if 0:
            print(row_df['biotypes'] + "\t" + row_df['biotype_of_div_file_path'] + "\t" + row_df['biotype_of_div_file'])
        row_df['biotype_of_div_file_handler'].write(  # add a new entry to the file
            my_species.name + "\t" + str(my_species.assembly) + "\t" + str(my_species.division) + "\t" + str(my_species.db) +
            "\t" + str(my_species.db_version) + "\t" + my_species.root_genes_path + "\t" + my_species.trunk_genes_path +
            "\t" + my_species.get_file_name_of_genes_of_species_and_biotype(str(g)) + "\t" +
            str(my_species.min_number_of_genes_per_biotype) + "\t" + str(number_of_genes) +
            "\t" + str(g) + "\n")

        # save entry in division_biotype_species file main output file
        div = ""
        if str(my_species.division)=="nan" and c.DIVISION=="vertebrates":
            div = c.DIVISION
        else:
            div = str(division)
        division_biotype_species_f_h.write(
            div + "\t" + str(g) + "\t" + str(my_species.name) + "\t" + str(my_species.assembly) + "\t" +
            str(my_species.db) + "\t" + str(my_species.db_version) + "\t" +
            my_species.root_annotation_path + "\t" + my_species.trunk_annotation_path + "\t" + my_species.annotation_file +
            "\t" + my_species.root_biotype_distribution_path  + "\t" + my_species.trunk_biotype_distribution_path + "\t" +
            c.LOG_BIOTYPE_DISTRIBUTION_SPECIES_FILE.replace(c.BIOTYPE_DISTR_OUTPUT_INPUT_PATH, "") + "\t" +
            my_species.root_genes_path + "\t" + my_species.trunk_genes_path + "\t" +
            my_species.get_file_name_of_genes_of_species_and_biotype(str(g)) + "\t" +
            str(my_species.min_number_of_genes_per_biotype) + "\t" + str(df_genes_of_biotype.shape[0]) + "\t" +
            str(df_all_genes.shape[0]) + "\t" + str(df_species.shape[0]) + "\n")
        sys.exit()
    if 1:
        total_time = time() - start_time
        print("...total time: %.10f seconds." % total_time)

# close opened files
f.close()  # the biotype_distribution_of_the_division file
for index, row in df_biotypes_of_div.iterrows(): # each biotype_of_division file
    row['biotype_of_div_file_handler'].close()
division_biotype_species_f_h.close()