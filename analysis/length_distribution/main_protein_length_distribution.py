from lib import constants as c
from lib_analysis import constants_analysis as ca
from lib import statistics as EM_stat
from lib import plot

import pandas as pd
import numpy as np
import os
from time import time
import sys

ENTITY = "protein"
BOOL_PLOT_OR_TABLE = "PLOT"   # if TABLE only create the table (no plot generation)
                               # if PLOT only generate the plots (png)

def main():
    start_time = time()

    # Protein file (stat_description) These are all the species (19,554)
    # 330 archaea, 7997 bacteria, 1588 eukaryota, 9939 viruses
    #
    prots_with_path = ca.OUTPUT_INPUT_FILES_PATH + ca.SOME_STATISTICS_PROTEINS_PATH_NAME + \
                      "/stat_description.protein.uniprot_reference_proteome.tsv"

    df_prots_stat_desc = pd.read_csv(prots_with_path, index_col=None, header=0, sep="\t")
    if 1: # protects from writting
        print(prots_with_path)
        pd.options.display.max_columns=None
        print(df_prots_stat_desc.shape)
        print(df_prots_stat_desc.head(3))
        print(df_prots_stat_desc.columns)
        print(str(len(pd.unique(df_prots_stat_desc['species']))) + " diff species")
        sys.exit()

    stat_col_names = ["species", "proteome_id", "tax_id", "superregnum", "num_prot_cod_genes", "uniprot_fasta_file",
                     "kurtosis", "skew", "stat_shapiro_wilk", "p_shapiro_wilk",
                     "stat_agostino_pearson", "p_agostino_pearson", "stat_kolmogorov_smirnov", "p_kolmogorov_smirnov"]
    df_genes_stat_logNorm = pd.DataFrame(columns=stat_col_names)

    #
    # iterate through each row (species)
    for i in range(len(df_prots_stat_desc)):
        annotation_name = df_prots_stat_desc.loc[i, "species"]
        superregnum = df_prots_stat_desc.loc[i, "superregnum"]
        if superregnum == "viruses": # I do not do viruses because they have not enough prots.
            continue

        list_id = [df_prots_stat_desc.loc[i, "species"], df_prots_stat_desc.loc[i, "proteome_id"],
                   df_prots_stat_desc.loc[i, "tax_id"], df_prots_stat_desc.loc[i, "superregnum"],
                   df_prots_stat_desc.loc[i, "num_prot_cod_genes"], df_prots_stat_desc.loc[i, "uniprot_fasta_file"]]
        if i % 20 == 0:
            partial_time = time() - start_time
            print(str(i) + "...partial time: %.10f seconds." % partial_time)
            print(list_id)

        protsLength_file_with_path = c.OUT_DATA_LOCAL_PATH_ROOT + \
                                     df_prots_stat_desc.loc[i, "uniprot_fasta_file"].replace(".fasta.gz", ".length.tsv")
        png_out_file_with_path=protsLength_file_with_path.replace(".tsv",".png")
        if 0:
            print(protsLength_file_with_path)
        if BOOL_PLOT_OR_TABLE == "PLOT":
            # the png file to save the plot of the distribution
            print(png_out_file_with_path)
            if os.path.isfile(png_out_file_with_path):
                continue

        # get all the proteins lengths (for the specie proteome)
        df_proteins = pd.read_csv(protsLength_file_with_path, sep="\t")
        df_proteins["log10ofLength"] = np.log10(df_proteins.length)
        if 0:
            pd.options.display.max_columns = None
            print(df_proteins.head(5))
            print(df_proteins.tail(5))
            sys.exit()

        stats_of_length = df_proteins[["length", "log10ofLength"]].describe()
        if 0:
            print(stats_of_length)
            print(type(stats_of_length["log10ofLength"]))
            print(stats_of_length["log10ofLength"]["count"].astype(int).item())
            sys.exit()
        if 0:
            print(stats_of_length["length"]["std"])
            sys.exit()

        data = df_proteins["log10ofLength"]
        #
        # is a log-norm distribution?
        kurtosis, skew,stat_shapiro_wilk, p_shapiro_wilk, stat_agostino_pearson, p_agostino_pearson, \
        stat_kolmogorov_smirnov, p_kolmogorov_smirnov = EM_stat.does_data_follow_a_normal_distr(data, 0)  # 1: print results

        text__is_a_logNorm = "Kurtosis:%.2f Skewness:%.2f\nShap-Wilk:\n %.3f p_val=%.3f\nAgost-Pears:\n %.3f p_val=%.3f\n" \
                             "Kol_Smirn:\n %.3f p_val=%.3f\n" % (kurtosis, skew, stat_shapiro_wilk, p_shapiro_wilk,
                             stat_agostino_pearson, p_agostino_pearson, stat_kolmogorov_smirnov, p_kolmogorov_smirnov)

        list_stat = [kurtosis, skew, stat_shapiro_wilk, p_shapiro_wilk, stat_agostino_pearson, p_agostino_pearson,
                     stat_kolmogorov_smirnov, p_kolmogorov_smirnov]
        df_proteins_new_row_stat_logNorm = pd.DataFrame(dict(zip(stat_col_names, list_id + list_stat)),
                                                     columns=stat_col_names, index=[0])  # using zip
        if 0:
            pd.options.display.max_columns=None
            print(df_proteins_new_row_stat_logNorm.shape)
            print(df_proteins_new_row_stat_logNorm)
            sys.exit()

        df_genes_stat_logNorm = df_genes_stat_logNorm.append(df_proteins_new_row_stat_logNorm, ignore_index=True)

        # distribution of lengths vs. pdf
        if BOOL_PLOT_OR_TABLE == "PLOT":
            plot.plot_withMatplotLib_dfColumn_histogram_Log10x(data,
                                                               stats_of_length["log10ofLength"]["count"], superregnum, annotation_name,
                                                               text__is_a_logNorm, ENTITY,
                                                               1, 4,
                                                               1, png_out_file_with_path)

    df_genes_stat_logNorm.sort_values(by=["kurtosis", "skew"], ascending=True, inplace=True)
    if 0:
        pd.options.display.max_columns=None
        print(df_genes_stat_logNorm.head(5))

    # save the results in a file:
    if BOOL_PLOT_OR_TABLE == "TABLE":
        # i.e /Volumes/Wes/results/geneLength/outputInputFiles/analysis/some_statistics/stat_description/all.logNorm_stat.description.ensembl.tsv
        all_logNorm_out_file = ca.OUTPUT_INPUT_FILES_PATH + ca.SOME_STATISTICS_PROTEINS_PATH_NAME + \
                               "/all.logNormStat.description.reference_proteome.tsv"
        print(all_logNorm_out_file)
        df_genes_stat_logNorm.to_csv(all_logNorm_out_file, sep="\t", index=False)
    #
    # display the computational time used
    if 1:
        print(str(df_genes_stat_logNorm.shape[0]) + " species has been processed (logNorm statistics)")
        total_time = time() - start_time
        print("...total time: %.10f seconds."%total_time)


if __name__=="__main__":
    main()