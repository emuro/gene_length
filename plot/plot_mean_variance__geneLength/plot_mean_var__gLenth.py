from lib_analysis import constants_analysis as c
import pandas as pd
import numpy as np
import math

import sys
#import plotnine as p9
from plotnine import *
from scipy import stats

BOOL_REVERSE = 0  # check I do not see much difference if I reverse the order of the divisions
DIVISIONS = ["viruses", "bacteria", "protist", "fungi", "plants", "metazoa", "vertebrates"]
if BOOL_REVERSE:
    DIVISIONS.reverse()
DBS = ["ensemblgenomes", "ensemblgenomes", "ensemblgenomes", "ensemblgenomes", "ensemblgenomes", "ensemblgenomes",
       "ensembl"] # for viruses, bateria, ...respectively
if BOOL_REVERSE:
    DBS.reverse()
ENSEMBL_VERSIONS = ["101", "49", "49", "49", "49", "49", "98"] # for viruses, bacteria, ...respectively
if BOOL_REVERSE:
    ENSEMBL_VERSIONS.reverse()

df = pd.DataFrame()
if 1:  # all genes of all species
    for i in range(len(DIVISIONS)):
        stat_file = c.OUTPUT_INPUT_FILES_PATH + c.SOME_STATISTICS_PATH_NAME + "/" + \
                        "stat_description." + DIVISIONS[i] + "." + DBS[i] + "." + str(ENSEMBL_VERSIONS[i]) + ".tsv"
        print(stat_file)
        df_aux = pd.read_csv(stat_file, sep="\t")
        df_aux["Division"] = DIVISIONS[i]
        print(df_aux[["mean","var"]].describe())

        if (i == 0):
            df = df_aux[:]
        else:
            df = df.append(df_aux, ignore_index=True)
        if 0:
            print(df.info())
            print(df[["mean", "var"]].describe())
else:  # tax_id merged (genes and prots). Species with well annotated genes and proteins.
    stat_file = "/Volumes/Wes/results/geneLength/outputInputFiles/analysis/some_statistics/stat_description/taxid_merged/stat_description.taxid_merged.ensembl_and_ref_proteome.tsv"
    print(stat_file)
    df = pd.read_csv(stat_file, sep="\t")
    df = df.rename(columns = {'genes_trunk_genes_path':'trunk_genes_path', 'genes_mean':'mean', 'genes_var':'var', 'merged_division_superregnum':'Division'} )


#
# PLOT
#######
df2plot = df.copy()
#df2plot = df2plot[df2plot.Division!="archaea"]  # possibility to filter manually some species
#df2plot = df2plot[df2plot.Division!="bacteria"]
df2plot = df2plot[df2plot.Division!="viruses"]
if 0:
    pd.options.display.max_columns=None
    print(df2plot.head(9))
    print(df2plot.shape[0])
    #my_list = df2plot["var"].to_list()
    #print(my_list)  # check that there are no nan value
    sys.exit()


BOOL_FILTER_GENES = False  # True
if BOOL_FILTER_GENES: # FILTER depending on the length distribution: Kurtosis, Skew
    # i.e /Volumes/Wes/results/geneLength/outputInputFiles/analysis/some_statistics/stat_description/all.logNorm_stat.description.ensembl.tsv
    all_logNorm_out_file = c.OUTPUT_INPUT_FILES_PATH + c.SOME_STATISTICS_PATH_NAME + \
                           "/all.logNormStat.description.ensembl.tsv"
    df_genes_stat_logNorm = pd.read_csv(all_logNorm_out_file, sep="\t")
    df_genes_stat_logNorm["new_kurtosis"] = df_genes_stat_logNorm["kurtosis"].abs()
    df_genes_stat_logNorm["new_skew"] = abs(df_genes_stat_logNorm["skew"])
    #
    # condition to filter based in how good fits the log-norm
    if 1:
        cond = df_genes_stat_logNorm["new_kurtosis"] >= 0.1
        df_genes_stat_logNorm.drop(df_genes_stat_logNorm[cond].index, inplace=True)
        cond2 = df_genes_stat_logNorm["new_skew"] >= 0.1
        df_genes_stat_logNorm.drop(df_genes_stat_logNorm[cond2].index, inplace=True)
    print(df_genes_stat_logNorm.shape)
    #
    if 1:
        print(all_logNorm_out_file)
        pd.options.display.max_columns=None
        print(df_genes_stat_logNorm.shape)
        print(df_genes_stat_logNorm.head(3))
        print(df_genes_stat_logNorm.columns)
        df_genes_stat_logNorm.to_csv("~/tmp/ks.tsv")
        print(str(len(pd.unique(df_genes_stat_logNorm['species'])))+" diff species")
        sys.exit()

    cond1 = df2plot["trunk_genes_path"].isin(df_genes_stat_logNorm["trunk_genes_path"])
    ##print(cond1)
    df2plot = df2plot[cond1]
    ## df2plot = df2plot.drop(df2plot[cond1].index, inplace=True)
    if 1:
        pd.options.display.max_columns=None
        print(df2plot.head(9))
        print(df2plot.shape[0])
        #sys.exit()

#calculate best fit line
slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(df2plot['mean']), np.log10(df2plot['var']))
txt = 'v = {:4.4} * m^{:4.4};   R^2= {:2.4f}'.format(10**intercept, slope, r_value**2)  # format text
p = (ggplot(df2plot, aes("mean", "var", color="Division")) + geom_point(size=0.05) +
     geom_smooth(method="lm", color="green", size=0.5, span=.8) +
     labs(title="Ensembl (protein coding gene length)")
     + scale_x_log10(breaks=[10 ** power for power in range(6)],
                   limits=[min(df2plot["mean"].to_list())/2, 2*max(df2plot["mean"].to_list())]
                   ) +
     scale_y_log10(breaks = [10**power for power in range(13)],
                   limits = [min(df2plot["var"].to_list())/2, 2*max(df2plot["var"].to_list())])
    ) +  theme(legend_position=(0.8,0.3), legend_key_size=5, legend_background=element_rect(fill='grey', alpha=0.01)) + \
    annotate('text', x=0.03*max(df2plot["mean"].to_list()), y=0.85*max(df2plot["var"].to_list()),
             label=txt, size=7, color="black")

print(p) # for some reason, I have to print manually the plot