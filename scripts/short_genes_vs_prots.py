from lib_analysis import constants_analysis as c
import sys
# Import the libraries
import matplotlib.pyplot as plt

import pandas as pd


def df_merged_filter_numOfProtsOrGenes(df_merged, num_threshold):
    """
         From the genes associated to proteines: ensembl-uniprot. Eliminate
         the species when the number of genes or proteins <= num_threshold

         Parameters
         ----------
         df_merged: df,
             the table with species (genes associated to proteins; ensembl to uniprot)
         num_threshold, int
             minimum number of genes or proteins
         Output:
         ------
         df_filtered,
             the table with the species filtered
    """
    BOOL_SHOW_COUNTS = 0
    if BOOL_SHOW_COUNTS:
        print("\t" + str(df_merged.shape[0]) + "\tbefore filtering")
    df_filtered = df_merged.loc[df_merged['prots_count']>num_threshold]
    if BOOL_SHOW_COUNTS:
        print("\t" + str(df_filtered.shape[0]) + "\tafter filtering prots")
    df_filtered = df_filtered.loc[df_filtered['genes_count']>num_threshold]
    if BOOL_SHOW_COUNTS:
        print("\t" + str(df_filtered.shape[0]) + "\tafter filtering genes")
    return  df_filtered


def df_merged_filter_diff_counts_genes_prots(df_merged,num_threshold_l,num_threshold_h):
    """
         From the genes associated to proteines: ensembl-uniprot. Eliminate
         the species when the number of genes or proteins:
            num_threshold_l <= X <= num_threshold_g

         Parameters
         ----------
         df_merged: df,
             the table with species (genes associated to proteins; ensembl to uniprot)
         num_threshold_l and num_threshold_h, int
             lower and higher thresholds
         Output:
         ------
         df_filtered,
             the table with the species filtered
    """
    BOOL_SHOW_COUNTS = 0
    if BOOL_SHOW_COUNTS:
        print("\t" + str(df_merged.shape[0]) + "\tbefore filtering")
    df_merged["diff_prots_genes"] = df_merged['prots_count']-df_merged['genes_count']
    df_merged["abs_diff"] = abs(df_merged['prots_count']-df_merged['genes_count'])
    df_merged["ratio_prots_genes"] = df_merged['prots_count']/df_merged['genes_count']
    if 0: # filt by the abs(diff)
        print( sorted(df_merged["abs_diff"].tolist()) )
        df_filtered = df_merged.loc[df_merged['abs_diff']<num_threshold]
        print(df_filtered.groupby("merged_division_superregnum").count().genes_species)
        print(df_filtered.groupby("merged_division_superregnum").count().genes_species.sum())
    if 1: # ratio
        df_filtered = df_merged.loc[df_merged['ratio_prots_genes']>num_threshold_l]
        df_filtered = df_filtered.loc[df_filtered['ratio_prots_genes']<num_threshold_h]
    return  df_filtered


def df_merged_filter_ratio_length_3xprots_genes(df_merged, num_threshold):
    """
         From the genes associated to proteines: ensembl-uniprot. Eliminate
         the species when the lenght of the genes >= 3*proteins

         Parameters
         ----------
         df_merged: df,
             the table with species (genes associated to proteins; ensembl to uniprot)
         num_threshold, int
             minimum number of genes or proteins
         Output:
         ------
         df_filtered,
             the table with the species filtered
    """
    BOOL_SHOW_COUNTS = 1
    if BOOL_SHOW_COUNTS:
        print("\t" + str(df_merged.shape[0]) + "\tbefore filtering")
    df_merged["ratio_mean_prots_genes"] = 3*df_merged['prots_mean']/df_merged['genes_mean']
    if 1: #ratio
        df_filtered = df_merged.loc[df_merged['ratio_mean_prots_genes']<num_threshold]
    return  df_filtered




df_all_merged =  pd.read_csv(c.LOG_SOME_STATISTICS_TAXID_MERGED_FILE, sep="\t")
df_aux = df_all_merged.copy()

if 1:
    print(df_aux.groupby("merged_division_superregnum").count().genes_species)
    print(df_aux.groupby("merged_division_superregnum").count().genes_species.sum())
    print()

if 1:
    df_filt_ratioMean = df_merged_filter_ratio_length_3xprots_genes(df_aux, 1.0)
    df_aux = df_filt_ratioMean.copy()
    print(df_aux.groupby("merged_division_superregnum").count().genes_species)
    print(df_aux.groupby("merged_division_superregnum").count().genes_species.sum())
    print()

if 1:
    df_filt_num = df_merged_filter_numOfProtsOrGenes(df_aux, 500)
    df_aux = df_filt_num.copy()
    print(df_aux.groupby("merged_division_superregnum").count().genes_species)
    print(df_aux.groupby("merged_division_superregnum").count().genes_species.sum())
    print()

if 1:
    df_filt_diffCount = df_merged_filter_diff_counts_genes_prots(df_aux, 0.95, 1.05)
    df_aux = df_filt_diffCount.copy()
    print(df_aux.groupby("merged_division_superregnum").count().genes_species)
    print(df_aux.groupby("merged_division_superregnum").count().genes_species.sum())
    print()


"""
if 0: # weird
    weird_path = "/Volumes/Wes/results/geneLength/outputInputFiles/analysis/some_statistics/stat_description/taxid_merged/problema_genesCortos_frente_longitudProteina/"
    #weird_file = weird_path + "anomalos_3prot_great_gene_G1.05.txt"      # 68
    weird_file = weird_path + "anomalos_3prot_great_gene_withHeader.txt"  # 989
else:
    all_file = c.LOG_SOME_STATISTICS_TAXID_MERGED_FILE
    weird_file = all_file

if 0:
    df_weird = pd.read_csv(weird_file, sep="\t")
    df_weird["ratio_weird"] = (3* df_weird.prots_mean)/df_weird.genes_mean
    print(df_weird)
    print(df_weird.columns)

    # matplotlib histogram
    plt.hist(df_weird['ratio_weird'], color = 'blue', edgecolor = 'black', bins = int(100)) #
    # Add labels
    plt.title('989 weird lengths: anomalos_3prot_great_gene_withHeader.txt')
    #plt.title('7675 merged')
    #plt.title('68 merged > 5%merged')
    #
    plt.ylim(0, 100)
    plt.xlim(1.01, 2.75)
    #
    plt.xlabel('3*prots_mean/genes_mean')
    plt.ylabel('num')
    #
    plt.show()
"""