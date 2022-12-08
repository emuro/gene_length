from lib_analysis import constants_analysis as c

import pandas as pd
import numpy as np
import math
from scipy import stats
import sys
from plotnine import *
from mizani.formatters import scientific_format

#Genes (Ensembl): taxid
#
df_merged_taxid = pd.read_csv(c.LOG_SOME_STATISTICS_TAXID_MERGED_FILE, sep="\t")
if 0:
    print(c.LOG_SOME_STATISTICS_TAXID_MERGED_FILE)
    pd.options.display.max_columns=None
    print(df_merged_taxid.head(3))
    print(df_merged_taxid.shape)
    print(str(len(pd.unique(df_merged_taxid['tax_id']))) + " diff tax_id")
    print(str(len(pd.unique(df_merged_taxid['genes_species']))) + " diff species (ensembl)")
    sys.exit()

df2plot = df_merged_taxid.copy()
#df2plot = df2plot[df2plot.merged_division_superregnum == "vertebrates"] # protist (my mistake not protists)

#
# mean vs. var (protein or gene)
if 0:
    col_x = "genes_mean"  # prots_mean
    col_y = "genes_var"   # prots_var
    x_lab = "mean (nt)"   # "mean (aa)"
    y_lab = "var"
    legends_by = "merged_division_superregnum"
    title = "Ensembl (protein coding gene length)" # "Uniprot (reference Proteomes)"

    #calculate best fit line
    slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(df2plot[col_x]),np.log10(df2plot[col_y]))
    #format text
    txt = 'v = {:4.4} * m^{:4.4};   R^2= {:2.4f}'.format(10**intercept, slope, r_value**2)
    p = (ggplot(df2plot, aes(col_x, col_y, color=legends_by)) + geom_point(size=0.75)+
         geom_smooth(method="lm", color="green", size=0.5, span=.8)+
         labs(title=title, x=x_lab, y=y_lab)
         + scale_x_log10(#breaks=[10 ** power for power in range(6)],
                         limits=[min(df2plot[col_x].to_list())/2, 2*max(df2plot[col_x].to_list())]
                         ) +
         scale_y_log10(#breaks = [10**power for power in range(13)],
                       limits = [min(df2plot[col_y].to_list())/2,2*max(df2plot[col_y].to_list())],
                       labels=scientific_format(digits=2))
         ) +  theme(legend_position=(0.8,0.3), legend_key_size=5, legend_background=element_rect(fill='grey', alpha=0.01)) + \
        annotate('text',x=0.2*max(df2plot[col_x].to_list()),y=0.85*max(df2plot[col_y].to_list()), # prots: 0.175,0.85
                 label=txt,size=7,color="black")

    print(p) # for some reason, I have to print my plot

#
# compare the mean of length (prots vs. genes)
if 0:
    df2plot["nt_prots_mean"] = 3 * df2plot["prots_mean"]
    df2plot["ratio_mean"] = df2plot["genes_mean"]/df2plot["nt_prots_mean"]
    #
    legends_by = "merged_division_superregnum"
    title = "Uniprot (ref. proteomes) vs Ensembl (prot. coding genes)"
    col_x="nt_prots_mean"  # genes_mean
    x_lab="proteins (nt)"
    if 1: # ratio in y_axis
        col_y = "ratio_mean"
        y_lab =  "ratio gene/prot (nt)"
    else: # pure length in y_axis
        col_y="genes_mean"
        y_lab="genes (nt)"

    p = (ggplot(df2plot, aes(col_x, col_y, color=legends_by)) + geom_point(size=0.75)+
         labs(title=title, x=x_lab, y=y_lab)
         )+  theme(legend_position=(0.375,0.7), legend_key_size=5, legend_background=element_rect(fill='grey', alpha=0.01))

    print(p) # for some reason, I have to print my plot

# Me pide Fer...
# compare the mean of length (prots vs. genes)
if 0:
    df2plot["nt_prots_mean"] = 3 * df2plot["prots_mean"]
    df2plot["ratio_mean"] = df2plot["genes_mean"]/df2plot["nt_prots_mean"]
    #
    legends_by = "merged_division_superregnum"
    title = "Uniprot (ref. proteomes) vs Ensembl (prot. coding genes)"
    col_x="genes_mean"  # genes_mean
    x_lab="mean (genes; nt)"
    # ratio in y_axis
    col_y = "ratio_mean"
    y_lab =  "ratio mean genes/prots (nt)"

    p = (ggplot(df2plot, aes(col_x, col_y, color=legends_by)) + geom_point(size=0.75)+
         labs(title=title, x=x_lab, y=y_lab) #+ xlim(0, 2500) + ylim(0, 2.5) #20000, 10
         +scale_x_log10(  # breaks=[10 ** power for power in range(6)],
                limits=[min(df2plot[col_x].to_list())/2,2*max(df2plot[col_x].to_list())]
            )+
         scale_y_log10(  # breaks = [10**power for power in range(13)],
             limits=[min(df2plot[col_y].to_list())/2,2*max(df2plot[col_y].to_list())],
             labels=scientific_format(digits=2))
         )+  theme(legend_position=(0.375,0.7), legend_key_size=5, legend_background=element_rect(fill='grey', alpha=0.01))

    print(p) # for some reason, I have to print my plot


#
# compare the var of length (prots vs. genes)
if 1:
    df2plot["nt_prots_var"] = df2plot["prots_var"]
    df2plot["ratio_var"] = df2plot["genes_var"]/df2plot["nt_prots_var"]
    #
    legends_by = "merged_division_superregnum"
    title = "Uniprot (ref. proteomes) vs Ensembl (prot. coding genes)"
    col_x="prots_var"  # genes_mean
    x_lab="proteins (var)"
    if 1: # ratio in y_axis
        col_y = "ratio_var"
        y_lab =  "ratio var genes/prots (nt)"
    else: # pure length in y_axis
        col_y="genes_var"
        y_lab="genes (var)"

    p = (ggplot(df2plot, aes(col_x, col_y, color=legends_by)) + geom_point(size=0.75)+
         labs(title=title, x=x_lab, y=y_lab)
         + scale_x_continuous(labels=scientific_format(digits=2))
         + scale_y_continuous(labels=scientific_format(digits=2))
         ) +  theme(legend_position=(0.7,0.7), legend_key_size=5, legend_background=element_rect(fill='grey', alpha=0.01))

    print(p) # for some reason, I have to print my plot