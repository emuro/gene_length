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
#
import pandas as pd
import sys
import glob
import re


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



# Genes (Ensembl): taxid
#
ensembl_taxid_file = c.OUTPUT_INPUT_FILES_PATH + c.ENSEMBL_TAXID_FILE_NAME
df_ensembl = pd.read_csv(ensembl_taxid_file, sep="\t")
if 0: # show the filtered ones in order
    print(df_ensembl.loc[df_ensembl['taxonomy_id'] == 8090])
    sys.exit()
if 0:
    print(ensembl_taxid_file)
    pd.options.display.max_columns=None
    print(df_ensembl.head(3))
    print(df_ensembl.shape)
    print( str(len( pd.unique(df_ensembl['taxonomy_id']))) + " diff tax_id")
    print( str(len(pd.unique(df_ensembl['species']))) + " diff species")
    sys.exit()
list_ensembl_taxid = df_ensembl['taxonomy_id'].to_list()
if 0:
    print(len(list_ensembl_taxid)) # redundancz of taxids: before uniq
list_ensembl_taxid = list(set(list_ensembl_taxid)) # avoid redundancy, there is
list_ensembl_taxid.sort(key=int)
if 0:
    print(list_ensembl_taxid)
    print(len(list_ensembl_taxid)) # redundancy of taxids: after uniq
    print(type(list_ensembl_taxid))
    sys.exit()


# Ensembl file (stat_description)
#
files_with_path = c.OUTPUT_INPUT_FILES_PATH + c.SOME_STATISTICS_PATH_NAME + "/stat_description.*.tsv"
genes_all_files = glob.glob(files_with_path)
if 0:
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
if 0:
    pd.options.display.max_columns=None
    print(df_genes_stat_desc.shape)
    print(df_genes_stat_desc.head(3))
    print(df_genes_stat_desc.columns)
    print(str(len(pd.unique(df_genes_stat_desc['species'])))+" diff species")
    sys.exit()

# Reference proteomes (Uniprot): taxid & prots.
#
stat_file = c.LOG_SOME_STATISTICS_REFERENCE_PROTEOME_FILE
df_prot_stat_desc = pd.read_csv(stat_file, sep="\t")
df_ref_prot = df_prot_stat_desc.copy()
if 0:
    print(stat_file)
    pd.options.display.max_columns=None
    print(df_prot_stat_desc.head(1))
    print(df_prot_stat_desc.shape)
    print(len( pd.unique(df_prot_stat_desc['tax_id']))) # checked: there is no redundancy
    sys.exit()
list_prot_taxid = df_prot_stat_desc['tax_id'].to_list()
list_prot_taxid = list(set(list_prot_taxid)) # avoid redundancy, there is
list_prot_taxid.sort(key=int)
if 0:
    print(list_prot_taxid)
    print(len(list_prot_taxid))
    sys.exit()

# Intersection of taxids: Ensembl and Uniprot
# Non redundant
list_intersect = []
for g in list_ensembl_taxid:
    if g in list_prot_taxid:
        list_intersect.append(g)
if 0:
    print(list_intersect)
    print(len(list_intersect))
    print(type(list_intersect))
    sys.exit()

# Intersection applied to Uniprot (taxid)
#
df_ref_prot_just_filtered = df_ref_prot[df_ref_prot['tax_id'].isin(list_intersect)]
df_ref_prot_filtered = df_ref_prot_just_filtered[["tax_id", "species", "proteome_id", "superregnum"]]
df_ref_prot_filtered = df_ref_prot_filtered.rename(columns={"species": "prot_species"})
df_ref_prot_filtered = df_ref_prot_filtered.sort_values("tax_id")
df_ref_prot_filtered.drop_duplicates(subset="tax_id", keep="first", inplace=True) # there was no redunt species here
if 0:
    pd.options.display.max_columns=None
    print(df_ref_prot_filtered.head(3))
    print(df_ref_prot_filtered.shape)
    print(len( pd.unique(df_ref_prot_filtered['tax_id']))) # checked: there is no redundancy
    sys.exit()
if 0:
    df_ref_prot_filtered.to_csv("/Users/enriquem.muro/tmp/prot.tsv", sep="\t") # for testing

# Intersection applied to ensembl (taxid)
#
df_ensembl_just_filtered = df_ensembl[df_ensembl['taxonomy_id'].isin(list_intersect)]  # maintains the order of df_ensembl
df_ensembl_just_filtered = df_ensembl_just_filtered.drop_duplicates(subset=["taxonomy_id"], keep="first") # the first species by taxonomy_id
if 0: # show the filtered ones in order
    print(df_ensembl_just_filtered.loc[df_ensembl['taxonomy_id'] == 8090])
    sys.exit()
df_ensembl_filtered = df_ensembl_just_filtered[["taxonomy_id", "species", "assembly", "division"]]
#df_ensembl_filtered = df_ensembl_filtered.sort_values("taxonomy_id") # NOOOO!!! todavia no!!! !!!!!!!!!
df_ensembl_filtered = df_ensembl_filtered.rename(columns={"taxonomy_id": "tax_id", "species": "genes_species"})
list_ensembl_filtered_species = list(set(df_ensembl_filtered['genes_species'].to_list()))
dict_ensembl_filtered_species = dict(zip(df_ensembl_filtered.genes_species, df_ensembl_filtered.tax_id))
if 0:
    pd.options.display.max_columns=None
    print(df_ensembl_filtered.head(2))
    print(df_ensembl_filtered.shape)
    print(len(list_ensembl_filtered_species))
    print(len( pd.unique(df_ensembl_filtered['tax_id'])))
    print(dict_ensembl_filtered_species)
    sys.exit()
if 0:
    df_ensembl_filtered.to_csv("/Users/enriquem.muro/tmp/genes.tsv", sep="\t") # for testing


# Intersection applied to ensembl (taxid)
# df_genes_stat_desc
df_genes_just_filtered = df_genes_stat_desc[df_genes_stat_desc['species'].isin(list_ensembl_filtered_species)]
df_genes_just_filtered = df_genes_just_filtered.reset_index(drop=True)
rule_a = (df_genes_just_filtered['species']=="caenorhabditis_elegans") & (df_genes_just_filtered['division']=='vertebrates')
rule_b = (df_genes_just_filtered['species']=="drosophila_melanogaster") & (df_genes_just_filtered['division']=='vertebrates')
rule_c = rule_a | rule_b
print(df_genes_just_filtered.shape)
for i, b in enumerate(rule_c):
    if b:
        df_genes_just_filtered.drop(i, axis=0, inplace=True)  # me ha costado bastante poder eliminar rows
df_genes_just_filtered = df_genes_just_filtered.reset_index(drop=True)
print(df_genes_just_filtered.shape)
df_genes_just_filtered['tax_id'] = df_genes_just_filtered['species'].map(dict_ensembl_filtered_species) # expectacular!
df_genes_just_filtered = df_genes_just_filtered.sort_values("tax_id")
df_genes_just_filtered.drop_duplicates(subset="tax_id", keep="first", inplace=True)  # eliminate redundant tax_id (but the first)
col = df_genes_just_filtered.columns # change the name of cols to genes_*, but tax_id
new_col = ["genes_" + c for c in col]
df_genes_just_filtered.columns = new_col
df_genes_just_filtered = df_genes_just_filtered.rename(columns={"genes_tax_id": "tax_id"})

# add the associated stat_description from the corresponding proteome
col = df_ref_prot_just_filtered.columns # change the name of cols to prots_*, but tax_id
new_col = ["prots_" + c for c in col]
df_ref_prot_just_filtered.columns = new_col
df_ref_prot_just_filtered = df_ref_prot_just_filtered.rename(columns={"prots_tax_id": "tax_id"})

# stat description merge: genes (Ensembl) vs. ref Proteome (Uniprot)
#
df_merged = pd.merge(df_genes_just_filtered, df_ref_prot_just_filtered, on='tax_id')
# improve the division/superregnum info (divide bacteria in Archaea and Bacteria, etc)
#
df_merged["merged_division_superregnum"] = df_merged["genes_division"]
cond1 = df_merged["prots_superregnum"] == "archaea" # there should be an easier way...but works fine!
cond2 = df_merged["merged_division_superregnum"] == "bacteria"
cond3 = cond1 & cond2
for i, b in enumerate(cond3):
    if b:
        df_merged.loc[i, "merged_division_superregnum"] = "archaea"
if 0:
    pd.options.display.max_columns=None
    print(df_merged.head(3))
    print(df_merged.shape)
    print(len(pd.unique(df_merged['tax_id']))) # checked: there is no redundancy
    sys.exit()


# FILTER to improve the annotation
df_aux = df_merged.copy()
if 1:
    print(df_aux.groupby("merged_division_superregnum").count().genes_species)
    print(df_aux.groupby("merged_division_superregnum").count().genes_species.sum())
    print()
#
if 1:
    df_filt_ratioMean = df_merged_filter_ratio_length_3xprots_genes(df_aux, 1.0)
    df_aux = df_filt_ratioMean.copy()
    print(df_aux.groupby("merged_division_superregnum").count().genes_species)
    print(df_aux.groupby("merged_division_superregnum").count().genes_species.sum())
    print()
#
if 1:
    df_filt_num = df_merged_filter_numOfProtsOrGenes(df_aux, 500)
    df_aux = df_filt_num.copy()
    print(df_aux.groupby("merged_division_superregnum").count().genes_species)
    print(df_aux.groupby("merged_division_superregnum").count().genes_species.sum())
    print()
#
if 1:
    df_filt_diffCount = df_merged_filter_diff_counts_genes_prots(df_aux, 0.95, 1.05)
    df_aux = df_filt_diffCount.copy()
    print(df_aux.groupby("merged_division_superregnum").count().genes_species)
    print(df_aux.groupby("merged_division_superregnum").count().genes_species.sum())
    print()


if 1:  # FALTA SALVAR RESULTADOS!
    df_aux.to_csv(c.LOG_SOME_STATISTICS_TAXID_MERGED_FILE, sep="\t", index=False)
