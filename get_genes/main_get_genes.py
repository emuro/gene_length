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
# and SKIP_FIRST_SPECIES (in this file)
#
# DIVISION = "bacteria"  # "viruses", "bacteria", "protist", "fungi", "plants", "metazoa", "vertebrates"
# There rest should adapt to this selection.
# Also take into account that the minimum number of genes, for protein_coding) is
# MIN_NUM_OF_GENES_PER_BIOTYPE = 0
# bacteria,
# OSError: [Errno 24] Too many open files: -> salinimonas_sp_hmf8227_gca_003143555
# Then SKIP_FIRST_SPECIES = 10234; 20468; 30702;
# SKIP_FIRST_SPECIES = 0 if not bacteria & starting from 1st bacteria
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


SKIP_FIRST_SPECIES = 0 # skip species: 10234; 20468; 30702; from the first species: 0

def get_filename_of_division_biotype_species(division="unknown_division",
                                             db="unknown_db", db_version="unknown_version"):
    """ i.e:  division_biotype_species.vertebrates.ensembl.98.tsv (nan is for vertebrates)
    """
    if str(division) == "nan" and c.DIVISION == "vertebrates":
        return "division_biotype_species." + c.DIVISION + "." + db + "." + str(db_version) + ".tsv"
    else:
        return "division_biotype_species." + str(division) + "." + db + "." + str(db_version) + ".tsv"


def get_header_for_division_biotype_species_out_file():
    """ get the header for the main output file: division_biotype_species file
    """
    return "division\tbiotype\tspecies\tassembly\tdb\tdb_version\troot_annotation_path\ttrunk_annotation_path\t" \
           "annotation_file\troot_biotype_distribution_path\ttrunk_biotype_distribution_path\t" \
           "biotype_distribution_file\troot_genes_path\ttrunk_genes_path\tgenes_file\tmin_number_of_genes_per_biotype" \
           "\tnumber_of_genes_biotype\tnumber_of_genes_total\tnumber_of_species_in_division\n"


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


def main():
    start_time = time()
    #
    # Get the species info.
    # NOTE: the species from these file has been previously indexed (see main_index_ensembl_gtfFiles.py)
    # and the indexed file was accessed to obtain the biotype distribution of each species
    # (see main_biotype_distribution_of_species andm main_index_ensembl_gtfFiles.py)
    df_species = files.get_df_from_tsv_file(c.LOG_BIOTYPE_DISTRIBUTION_SPECIES_FILE) # species of the division
    for i in range(len(df_species)):  # Fix typo: iterate through each species and replace "//" by "/"
        df_species.loc[i, "trunk_biotype_distribution_path"] = \
            df_species.loc[i, "trunk_biotype_distribution_path"].replace("//", "/")
    if 0:
        pd.options.display.max_columns = None
        print(df_species.head(5))
        sys.exit("...jut got df_species with the biotype distribution")
    if 0:
        print(df_species.columns)

    #
    # prepare the main output file (division_biotype_species)
    division_biotype_species_out_file = get_filename_of_division_biotype_species(df_species.loc[i, "division"],
                                                                                 df_species.loc[i, "db"],
                                                                                 df_species.loc[i, "db_version"])
    pathlib.Path(c.GENES_OUTPUT_INPUT_PATH).mkdir(parents=True, exist_ok=True)  # just if it was not previously created
    division_biotype_species_f_h = open(c.GENES_OUTPUT_INPUT_PATH + division_biotype_species_out_file, "a")  # open file
    if not SKIP_FIRST_SPECIES:
        division_biotype_species_f_h.write(get_header_for_division_biotype_species_out_file())
    if 1:
        print(c.GENES_OUTPUT_INPUT_PATH + "\t" + division_biotype_species_out_file)


    for i in range(len(df_species)):  # iterate through each species
        if i < SKIP_FIRST_SPECIES:     # 10234: correcto -1 // ulimit -n 10240 (*4 vale)
            continue
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

            # save entry in division_biotype_species file main output file
            div = ""
            if str(my_species.division) == "nan" and c.DIVISION == "vertebrates":
                div = c.DIVISION
            else:
                div = str(my_species.division)
            division_biotype_species_f_h.write(
                div + "\t" + str(g) + "\t" + str(my_species.name) + "\t" + str(my_species.assembly) + "\t" +
                str(my_species.db) + "\t" + str(my_species.db_version) + "\t" +
                my_species.root_annotation_path + "\t" + my_species.trunk_annotation_path + "\t" +
                my_species.annotation_file +
                "\t" + my_species.root_biotype_distribution_path  + "\t" + my_species.trunk_biotype_distribution_path +
                "\t" + c.LOG_BIOTYPE_DISTRIBUTION_SPECIES_FILE.replace(c.BIOTYPE_DISTR_OUTPUT_INPUT_PATH, "") + "\t" +
                my_species.root_genes_path + "\t" + my_species.trunk_genes_path + "\t" +
                my_species.get_file_name_of_genes_of_species_and_biotype(str(g)) + "\t" +
                str(my_species.min_number_of_genes_per_biotype) + "\t" + str(df_genes_of_biotype.shape[0]) + "\t" +
                str(df_all_genes.shape[0]) + "\t" + str(df_species.shape[0]) + "\n")
        if 1:
            total_time = time() - start_time
            print("...total time: %.10f seconds." % total_time)

    # close opened files
    division_biotype_species_f_h.close()


if __name__ == "__main__":
    main()