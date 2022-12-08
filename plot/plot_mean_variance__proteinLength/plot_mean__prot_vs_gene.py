from lib_analysis import constants_analysis as c
import pandas as pd
import sys
import numpy as np
import glob

#import plotnine as p9
from plotnine import *
if 0:
    print(ggplot.__doc__)

from scipy import stats



df = pd.DataFrame()

# proteins
#
stat_file = c.LOG_SOME_STATISTICS_PROTEIN_FILE
df_prot_stat_desc = pd.read_csv(stat_file, sep="\t")
if 1:
    print(stat_file)
    pd.options.display.max_columns=None
    print(df_prot_stat_desc.head(9))
    print(df_prot_stat_desc.species.to_list())
    sys.exit()


# genes
#
files_with_path = c.OUTPUT_INPUT_FILES_PATH + c.SOME_STATISTICS_PATH_NAME + "/*.tsv"
genes_all_files = glob.glob(files_with_path)
if 0:
    print(genes_all_files)
    sys.exit()
li = []
for gene_file_name in genes_all_files:
    df = pd.read_csv(gene_file_name, index_col=None, header=0, sep="\t")
    li.append(df)
df_genes = pd.concat(li, axis=0, ignore_index=True)
df_genes = df_genes[["species", "count", "mean", "var" ]]

df_genes = df_genes[df_genes['species'].isin( df_prot_stat_desc.species.to_list() )]
df_genes = df_genes.drop_duplicates(subset="species", keep="last")
df_genes = df_genes.rename(columns={"count": "genes_count", "mean": "genes_mean", "var": "genes_var"})
df_genes = df_genes.reset_index(drop=True)
df_mix = pd.merge(df_prot_stat_desc, df_genes, how='inner')
df_mix = df_mix[["species", "count", "mean", "var", "genes_count", "genes_mean", "genes_var" ]]
df_mix = df_mix.rename(columns={"count": "proteins_count", "mean": "proteins_mean", "var": "proteins_var"})
if 0:
    pd.options.display.max_columns=None
    #rint(df_genes.shape)
    #print(df_prot_stat_desc.head(9))
    #print(df_genes.head(9))
    print(df_mix.head(9))
    sys.exit()


df2plot = df_mix.copy()

#calculate best fit line
slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(df2plot['proteins_mean']), np.log10(df2plot['genes_mean']))
#format text
txt = 'v = {:4.4} * m^{:4.4};   R^2= {:2.4f}'.format(10**intercept, slope, r_value**2)
print(txt)
p = (ggplot(df2plot, aes("proteins_mean", "genes_mean", color="species")) + geom_point(size=0.75) +
     geom_smooth(method="lm", color="green", size=0.2, span=.8)
     + labs(title="Proteins vs genes (comparison of means; model organisms)")
    ) +  theme(legend_position=(0.3,0.75), legend_key_size=5, legend_background=element_rect(fill='grey', alpha=0.01)) + \
    annotate('text', x=0.9*max(df2plot["proteins_mean"].to_list()), y=0.95*max(df2plot["genes_mean"].to_list()),
             label=txt, size=7, color="black")
if 1:
    print(p)


#calculate best fit line
slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(df2plot['proteins_var']), np.log10(df2plot['genes_var']))
#format text
txt = 'v = {:4.4} * m^{:4.4};   R^2= {:2.4f}'.format(10**intercept, slope, r_value**2)
print(txt)
p1 = (ggplot(df2plot, aes("proteins_var", "genes_var", color="species")) + geom_point(size=0.75) +
     geom_smooth(method="lm", color="green", size=0.2, span=.8)
     + labs(title="Proteins vs genes (comparison of variances; model organisms)")
    ) +  theme(legend_position=(0.3,0.75), legend_key_size=5, legend_background=element_rect(fill='grey', alpha=0.01)) + \
    annotate('text', x=0.9*max(df2plot["proteins_var"].to_list()), y=0.95*max(df2plot["genes_var"].to_list()),
             label=txt, size=7, color="black")
if 1:
    print(p1)

