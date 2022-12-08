# python3
# ################################################################## #
# main_some_statistics_proteins.py (C) March-April-2021 Mainz.
# Author: Enrique M. Muro
# ################################################################## #
#
# --------------------------------------------------------------------
# Project: geneLength
#
# Purpose:
#
# Important:
# Change:
# ################################################################## #
from lib_analysis import constants_analysis as c
from lib import EM_biopython_extras as bioEx
#
import pandas as pd
import gzip
from Bio import SeqIO
import numpy as np
from time import time
import os.path
import pathlib
import sys

def get_prot_length_description(df_prot,
                                species="", proteome_id="", tax_id="", superregnum="", num_prot_cod_genes="",
                                uniprot_fasta_file=""):
    """
        From a *.tsv file with headers, get a list from an specified column

        Parameters
        ----------
        df_prot: df
            length of all the proteins of a species
        species: str
            species name
        proteome_id: str
            proteome_ID of the species, as given by Uniprot
        tax_id: str
            taxonomical_ID of the species
        superregnum: str
            Bacteria, Viruses, Archaea or Eukaryota
        num_prot_cod_genes: int
            the number of prot cod. genes

        Output:
        ------
        df_description: df
            statistical description of the protein lengths for a species
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
    new_row = [species,  proteome_id, tax_id, superregnum, num_prot_cod_genes, uniprot_fasta_file] + list_n + list_log10 + list_log
    col_names = ["species", "proteome_id", "tax_id", "superregnum", "num_prot_cod_genes", "uniprot_fasta_file",
        "count", "mean", "std", "var", "min", "25perc", "50perc", "75perc", "max",
        "log10_mean", "log10_std", "log10_var", "log10_min", "log10_25perc", "log10_50perc", "log10_75perc", "log10_max",
        "log_mean", "log_std", "log_var", "log_min", "log_25perc", "log_50perc", "log_75perc", "log_max"]

    df_genes_description = pd.DataFrame(dict(zip(col_names, new_row)), columns=col_names, index=[0])  # using zip
    return df_genes_description

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
    # Get the species info and the biotype distribution for each species
    ref_proteomes_file = c.LOG_SOME_STATISTICS_REFERENCE_PROTEOME_FILE  # c.OUTPUT_INPUT_FILES_PATH + c.REFERENCE_PROTEOMES_PATH_NAME + c.REFERENCE_PROTEOMES_FILE
    ### ref_proteomes_file = c.OUTPUT_INPUT_FILES_PATH + c.REFERENCE_PROTEOMES_PATH_NAME + c.REFERENCE_PROTEOMES_FILE
    df_ref_proteomes = pd.read_csv(ref_proteomes_file, sep="\t")
    if 1: # 1 to protect from saving/writting files
        print(ref_proteomes_file)
        pd.options.display.max_columns=None
        print(df_ref_proteomes.head(5))
        print(df_ref_proteomes.shape[0])
        sys.exit()


    col_names=["db","UniprotID","EntryName","ProteinName"]+["OS","OX","GN","PE","SV"]+["length"]
    for i in range(len(df_ref_proteomes)):  # iterate through each species
        if i % 50 == 0:
            print(str(i)  + "\t" + df_ref_proteomes.loc[i, "species"] + "\t" + df_ref_proteomes.loc[i, "proteome_id"])
            if 0:
                print("proteome_id" + "_" + "tax_id" + "_" + "superregnum" + "_" + "num_prot_cod_genes")
                print(df_ref_proteomes.loc[i, "proteome_id"] + "_" + str(df_ref_proteomes.loc[i, "tax_id"]) +
                      "_" + df_ref_proteomes.loc[i, "superregnum"] + "_" + str(df_ref_proteomes.loc[i, "num_prot_cod_genes"]))

        uniprot_file_with_path = c.BASE_PATH + c.DATA_PATH_NAME
        uniprot_path_ext = df_ref_proteomes.loc[i, "uniprot_fasta_file"]
        uniprot_file_with_path += uniprot_path_ext

        if 0:
            print(uniprot_file_with_path)
            print(uniprot_path_ext)
            sys.exit()

        df_prot_len=pd.DataFrame([], columns=col_names)
        new_dict_row={"db":[],"UniprotID":[],"EntryName":[],"ProteinName":[],
                      "OS":[],"OX":[],"GN":[],"PE":[],"SV":[],"length":[]}
        with gzip.open(uniprot_file_with_path, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                db, UniprotID, EntryName, ProteinName, OS, OX, GN, PE, SV = bioEx.parse_uniprot_header(record)
                # La forma mas eficiente que he encontrado
                new_dict_row["db"].append(db)
                new_dict_row["UniprotID"].append(UniprotID)
                new_dict_row["EntryName"].append(EntryName)
                new_dict_row["ProteinName"].append(ProteinName)
                new_dict_row["OS"].append(OS)
                new_dict_row["OX"].append(OX)
                new_dict_row["GN"].append(GN)
                new_dict_row["PE"].append(PE)#
                new_dict_row["SV"].append(SV)
                new_dict_row["length"].append(str(len(record.seq)))
                if 0:
                    print("\ndb:(%s)\tUniprotID:(%s)\tEntryName:(%s)\tProteinName:(%s)"%(db,UniprotID,EntryName,ProteinName))
                    print("OS:(%s)\tOX:(%s)\tGN:(%s)\tPE:(%s)\tSV:(%s)"%(OS,OX,GN,PE,SV))
                    sys.exit()
            #add rows: fastest way I found
            df_prot_len = pd.DataFrame.from_dict(new_dict_row, orient='index').transpose()
            df_prot_len['length']=df_prot_len['length'].astype(int)
            df_prot_len = df_prot_len.sort_values(by=["length"],)
            calculate_length_diff_to_closest_gene(df_prot_len)
            if 0:
                pd.options.display.max_columns=None
                print(df_prot_len.head(9))
                df_prot_len.to_csv("/Users/enriquem.muro/tmp/kk.tsv", sep="\t")
                sys.exit()
            # save_file
            out_fullPath = c.OUT_DATA_LOCAL_PATH_ROOT + os.path.dirname(uniprot_path_ext)
            out_filename = os.path.basename(uniprot_path_ext)
            out_filename = out_filename.replace(".fasta.gz", ".length.tsv")
            out_fullPath_filename = out_fullPath + "/" + out_filename
            if 0:
                print(out_fullPath)
                print(out_fullPath_filename)
                sys.exit()

            pathlib.Path(out_fullPath).mkdir(parents=True, exist_ok=True)
            df_prot_len.to_csv(out_fullPath_filename, sep="\t", index=False)

    if 1:
        print(str(df_ref_proteomes.shape[0]) + " proteomes has been processed (get_proteins, lenght, etc)")
        total_time = time() - start_time
        print("...total time: %.10f seconds." % total_time)


if __name__ == "__main__":
    main()
