from lib_analysis import constants_analysis as ca

import pandas as pd
import sys
import numpy as np
from plotnine import *
from scipy import stats


df = pd.DataFrame()

if 0:
    prot_stat_file = ca.LOG_SOME_STATISTICS_REFERENCE_PROTEOME_FILE
    df_prot_stat_desc = pd.read_csv(prot_stat_file, sep="\t")
    if 0:
        print(prot_stat_file)
        pd.options.display.max_columns=None
        print(df_prot_stat_desc.head(9))
        print(df_prot_stat_desc.shape)
        sys.exit()

prot_logNorm_stat_file = all_logNorm_out_file = ca.OUTPUT_INPUT_FILES_PATH + ca.SOME_STATISTICS_PROTEINS_PATH_NAME+ \
                         "/all.logNormStat.description.reference_proteome.tsv"
df_prot_logNorm_stat_desc = pd.read_csv(prot_logNorm_stat_file, sep="\t")
df_prot_logNorm_stat_desc = df_prot_logNorm_stat_desc[["species", "proteome_id", "tax_id", "superregnum", "kurtosis", "skew"]]
df_prot_logNorm_stat_desc = df_prot_logNorm_stat_desc.rename(columns = {'kurtosis':'proteins_kurtosis',
                                                                        'skew':'proteins_skewness'})
if 0:
    print(prot_logNorm_stat_file)
    pd.options.display.max_columns=None
    print(df_prot_logNorm_stat_desc.head(9))
    print(df_prot_logNorm_stat_desc.shape)
    sys.exit()

genes_logNorm_stat_file = all_logNorm_out_file = ca.OUTPUT_INPUT_FILES_PATH + ca.SOME_STATISTICS_PATH_NAME+ \
                         "/all.logNormStat.description.ensembl.tsv"
df_genes_logNorm_stat_desc = pd.read_csv(genes_logNorm_stat_file, sep="\t")
#df_genes_logNorm_stat_desc = df_genes_logNorm_stat_desc[["species", "assembly", "division", "trunk_genes_path", "kurtosis", "skew"]]
df_genes_logNorm_stat_desc = df_genes_logNorm_stat_desc[["trunk_genes_path", "kurtosis", "skew"]]
df_genes_logNorm_stat_desc = df_genes_logNorm_stat_desc.rename(columns = {'trunk_genes_path':'genes_trunk_genes_path',
                                                                          'kurtosis':'genes_kurtosis', 'skew':'genes_skewness'})
if 0:
    print(genes_logNorm_stat_file)
    pd.options.display.max_columns=None
    print(df_genes_logNorm_stat_desc.head(9))
    print(df_genes_logNorm_stat_desc.shape)
    sys.exit()


merged_stat_file = ca.LOG_SOME_STATISTICS_TAXID_MERGED_FILE
df_merged_stat_desc = pd.read_csv(merged_stat_file, sep="\t")
df_merged_stat_desc = df_merged_stat_desc[["genes_species", "genes_trunk_genes_path", "tax_id", "merged_division_superregnum",
                                           "prots_proteome_id", "genes_count", "prots_num_prot_cod_genes"]]
if 0:
    print(merged_stat_file)
    pd.options.display.max_columns=None
    print(df_merged_stat_desc.head(9))
    print(df_merged_stat_desc.shape)
    sys.exit()


# filter
cond1 = df_genes_logNorm_stat_desc["genes_trunk_genes_path"].isin(df_merged_stat_desc["genes_trunk_genes_path"])
df_genes_logNorm_stat_desc = df_genes_logNorm_stat_desc[cond1]
df_merged = pd.merge(df_merged_stat_desc, df_genes_logNorm_stat_desc, on='genes_trunk_genes_path')
#
cond2 = df_prot_logNorm_stat_desc["tax_id"].isin(df_merged_stat_desc["tax_id"])
df_prot_logNorm_stat_desc = df_prot_logNorm_stat_desc[cond2]
df_merged = pd.merge(df_merged, df_prot_logNorm_stat_desc, on='tax_id')
if 0:
    pd.options.display.max_columns=None
    print(df_merged.head(9))
    print(df_merged.shape)
    sys.exit()


df2plot = df_merged.copy()
df2plot = df2plot[["genes_species", "tax_id", "merged_division_superregnum", "genes_kurtosis", "genes_skewness",
                   "proteins_kurtosis", "proteins_skewness"]]
#

X = "plants"  # bacteria, archaea, fungi, protist, plants, metazoa, vertebrates
cond_X = df2plot["merged_division_superregnum"] == "metazoa"
df2plot = df2plot[cond_X]


#
cond2 = df_prot_logNorm_stat_desc["tax_id"].isin(df_merged_stat_desc["tax_id"])
df_prot_logNorm_stat_desc = df_prot_logNorm_stat_desc[cond2]
#df2plot = df2plot[df2plot.superregnum=="eukaryota"]
if 0:
    pd.options.display.max_columns=None
    print(df2plot.head(9))
    sys.exit()
print(df2plot.shape)

p = (ggplot(df2plot, aes("genes_kurtosis", "proteins_kurtosis", color="merged_division_superregnum")) + geom_point(size=0.5)
     + labs(title="Kurtosis")
     + geom_hline(yintercept=0)  # add one horizonal line
     + geom_vline(xintercept=0)  # add one vertical line
     #+ xlim(-1, 1) + ylim(-1,1)
     + xlim(-2,2)+ylim(-2,2)
    ) +  theme(legend_position=(0.3, 0.8), legend_key_size=5, legend_background=element_rect(fill='grey', alpha=0.01))

if 1:
    print(p)

p2 = (ggplot(df2plot, aes("genes_skewness", "proteins_skewness", color="merged_division_superregnum")) + geom_point(size=0.5)
     + labs(title="Skewness")
     + geom_hline(yintercept=0)  # add one horizonal line
     + geom_vline(xintercept=0)  # add one vertical line
     #+ xlim(-1,1) + ylim(-1,1)
     + xlim(-2,2)+ylim(-2,2)
    ) +  theme(legend_position=(0.3, 0.8), legend_key_size=2, legend_background=element_rect(fill='grey', alpha=0.01))
if 1:
    print(p2)
