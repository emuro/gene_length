from lib import constants as c
from lib_analysis import constants_analysis as ca
from lib import statistics as EM_stat
from lib import plot

# import constants as c
import pandas as pd
import glob
import re
import numpy as np
from time import time
import os
import sys

ENTITY = "gene"
BOOL_PLOT_OR_TABLE = "PLOT"   # if TABLE only create the table (no plot generation)
                              # if PLOT only generate the plots (png)

def main():
    start_time = time()

    # Ensembl file (stat_description) These are all the species.
    # bacteria:31943 fungi:1014 protist:237 vertebrates:227 metazoa:115 plants:96 viruses:1
    # 33633 and from those 32948 diff. species
    #
    files_with_path = ca.OUTPUT_INPUT_FILES_PATH + ca.SOME_STATISTICS_PATH_NAME + "/stat_description.*.tsv"
    genes_all_files = glob.glob(files_with_path)
    if 0:
        print(files_with_path)
        print(genes_all_files)
        sys.exit()
    li = []
    for gene_file_name in genes_all_files:
        division = re.search(r'stat_description\.(.*?).ensembl', gene_file_name).group(1)
        if 0:
            print(division)
        df = pd.read_csv(gene_file_name, index_col=None, header=0, sep="\t")
        df["division"] = division
        li.append(df)
    df_genes_stat_desc = pd.concat(li, axis=0, ignore_index=True)
    del df_genes_stat_desc['Unnamed: 0']
    if 1:  # 1 protect from running the whole code
        print(files_with_path)
        print(df_genes_stat_desc["division"].value_counts())
        pd.options.display.max_columns=None
        print(df_genes_stat_desc.shape)
        print(df_genes_stat_desc.head(3))
        print(df_genes_stat_desc.columns)
        df_genes_stat_desc__uniq = pd.unique(df_genes_stat_desc['species'])
        print(str(len(df_genes_stat_desc__uniq))+" diff species")
        #print(str(len(pd.unique(df_genes_stat_desc['species'])))+" diff species")
        df_sorted=df_genes_stat_desc.sort_values(['species'],ascending=[True])
        df_first=df_sorted.groupby('species').first().reset_index()
        print(df_first["division"].value_counts())
        sys.exit()


    stat_col_names = ["species", "assembly", "trunk_genes_path", "genes_file", "division",
                     "kurtosis", "skew", "stat_shapiro_wilk", "p_shapiro_wilk", "stat_agostino_pearson", "p_agostino_pearson",
                     "stat_kolmogorov_smirnov", "p_kolmogorov_smirnov"]
    df_genes_stat_logNorm = pd.DataFrame(columns=stat_col_names)

    #
    # iterate through each row (species)
    for i in range(len(df_genes_stat_desc)):
        annotation_name = df_genes_stat_desc.loc[i, "species"]
        reference_name = df_genes_stat_desc.loc[i, "assembly"]
        division = df_genes_stat_desc.loc[i, "division"]
        if c.BOOL_CHECK_SOME_SPECIES: #for the sake of a unique species analysis.
            if annotation_name.capitalize() not in c.ANNOTATION_NAMES_TO_CHECK:
                continue

        list_id = [df_genes_stat_desc.loc[i, "species"], df_genes_stat_desc.loc[i, "assembly"],
                   df_genes_stat_desc.loc[i, "trunk_genes_path"], df_genes_stat_desc.loc[i, "genes_file"],
                   df_genes_stat_desc.loc[i, "division"]]
        if i % 100 == 0:
            partial_time=time()-start_time
            print(str(i) + "...partial time: %.10f seconds." % partial_time)
            print(list_id)

        if BOOL_PLOT_OR_TABLE == "PLOT":
            # the png file to save the plot of the distribution
            png_out_path = c.OUT_DATA_LOCAL_PATH_ROOT+df_genes_stat_desc.loc[i,"trunk_genes_path"]
            png_out_file = df_genes_stat_desc.loc[i,"genes_file"].replace(".genes.",".geneLength_distrib.")
            png_out_file_with_path = png_out_path+png_out_file.replace(".tsv",".png")
            if os.path.isfile(png_out_file_with_path):
                continue

        gene_distr_file = c.OUT_DATA_LOCAL_PATH_ROOT + df_genes_stat_desc.loc[i,"trunk_genes_path"] + \
                          df_genes_stat_desc.loc[i, "genes_file"]
        if 0:
            print(gene_distr_file)

        # get all the gene lengths (for specie, gene_type)
        df_genes = pd.read_csv(gene_distr_file, sep="\t")
        df_genes["log10ofLength"] = np.log10(df_genes.length)
        if 0:
            pd.options.display.max_columns = None
            print(df_genes.head(5))
            print(df_genes.tail(5))
            sys.exit()

        stats_of_length = df_genes[["length", "log10ofLength"]].describe()
        if 0:
            print(stats_of_length)
            print(type(stats_of_length["log10ofLength"]))
            a = stats_of_length["log10ofLength"]["count"].astype(int).item()
            print(a)
            sys.exit()
        if 0:
            print(stats_of_length["length"]["std"])
            sys.exit()

        data = df_genes["log10ofLength"]
        #
        # is a log-norm distribution?
        kurtosis, skew,stat_shapiro_wilk, p_shapiro_wilk, stat_agostino_pearson, p_agostino_pearson, \
        stat_kolmogorov_smirnov, p_kolmogorov_smirnov = EM_stat.does_data_follow_a_normal_distr(data, 0)  # 1: print results

        text__is_a_logNorm = "Kurtosis:%.2f Skewness:%.2f\nShap-Wilk:\n %.3f p_val=%.3f\nAgost-Pears:\n %.3f p_val=%.3f\n" \
                             "Kol_Smirn:\n %.3f p_val=%.3f\n" % (kurtosis, skew, stat_shapiro_wilk, p_shapiro_wilk,
                             stat_agostino_pearson, p_agostino_pearson, stat_kolmogorov_smirnov, p_kolmogorov_smirnov)

        list_stat = [kurtosis, skew, stat_shapiro_wilk, p_shapiro_wilk, stat_agostino_pearson, p_agostino_pearson,
                     stat_kolmogorov_smirnov, p_kolmogorov_smirnov]
        df_genes_new_row_stat_logNorm = pd.DataFrame(dict(zip(stat_col_names, list_id + list_stat)),
                                                     columns=stat_col_names, index=[0])  # using zip
        if 0:
            pd.options.display.max_columns=None
            print(df_genes_new_row_stat_logNorm.shape)
            print(df_genes_new_row_stat_logNorm)
            sys.exit()

        df_genes_stat_logNorm = df_genes_stat_logNorm.append(df_genes_new_row_stat_logNorm, ignore_index=True)

        # distribution of lengths vs. pdf
        if BOOL_PLOT_OR_TABLE == "PLOT":
            # the png file to save the plot of the distribution
            png_out_path = c.OUT_DATA_LOCAL_PATH_ROOT + df_genes_stat_desc.loc[i, "trunk_genes_path"]
            png_out_file = df_genes_stat_desc.loc[i,"genes_file"].replace(".genes.", ".geneLength_distrib.")
            png_out_file_with_path = png_out_path + png_out_file.replace(".tsv", ".png")

            plot.plot_withMatplotLib_dfColumn_histogram_Log10x(data,
                                                               stats_of_length["log10ofLength"]["count"], division, annotation_name,
                                                               text__is_a_logNorm, ENTITY,
                                                               1, 7,
                                                               1, png_out_file_with_path)

        if 0: # THIS plot.plot_dfColumn_histogram_Log10x IS NOT READY
            plot.plot_dfColumn_histogram_Log10x(df=df_genes, col="log10ofLength", width_of_bin=0.1,
                                                count=int(stats_of_length["log10ofLength"]["count"]),
                                                mean=stats_of_length["log10ofLength"]["mean"],
                                                std=stats_of_length["log10ofLength"]["std"],
                                                min=stats_of_length["log10ofLength"]["min"],
                                                max=stats_of_length["log10ofLength"]["max"],
                                                plot_title=annotation_name, x_axis_label="length (log10 scaled)",
                                                bool_plot=0, plot_file_withPath="") #  plot: 0, 1
            sys.exit()

    df_genes_stat_logNorm.sort_values(by=["stat_kolmogorov_smirnov", "stat_shapiro_wilk", "stat_agostino_pearson"],
                                      ascending=False, inplace=True)
    if 0:
        pd.options.display.max_columns=None
        print(df_genes_stat_logNorm.head(5))

    # save the results in a file:
    if BOOL_PLOT_OR_TABLE == "TABLE":
        # i.e /Volumes/Wes/results/geneLength/outputInputFiles/analysis/some_statistics/stat_description/all.logNorm_stat.description.ensembl.tsv
        all_logNorm_out_file = ca.OUTPUT_INPUT_FILES_PATH + ca.SOME_STATISTICS_PATH_NAME + \
                               "/all.logNormStat.description.ensembl.tsv"
        print(all_logNorm_out_file)
        df_genes_stat_logNorm.to_csv(all_logNorm_out_file, sep="\t", index=False)
    #
    # display the computational time used
    if 1:
        print(str(df_genes_stat_logNorm.shape[0])+" species has been processed (logNorm statistics)")
        total_time=time()-start_time
        print("...total time: %.10f seconds."%total_time)


if __name__=="__main__":
    main()
