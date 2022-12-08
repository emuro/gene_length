from lib_analysis import constants_analysis as c
import pandas as pd
import sys
import numpy as np

#import plotnine as p9
from plotnine import *
##print(ggplot.__doc__)
from scipy import stats


df = pd.DataFrame()

stat_file = c.LOG_SOME_STATISTICS_REFERENCE_PROTEOME_FILE
if 1:
    print(stat_file)

df_prot_stat_desc = pd.read_csv(stat_file, sep="\t")
df_prot_stat_desc.dropna(inplace=True)
df_prot_stat_desc["Human_virus"]=False
if 0:
    pd.options.display.max_columns=None
    print(df_prot_stat_desc.head(9))
    sys.exit()

df_viruses = pd.read_csv("/Users/enriquem.muro/Desktop/viruses/human257_sorted_spNames_viruses__stat_description.protein.uniprot_reference_proteome.tsv", sep="\t", header=None)
#df_viruses = pd.read_csv("/Users/enriquem.muro/Desktop/viruses/human123_human257_sorted_spNames_viruses__stat_description.protein.uniprot_reference_proteome.tsv", sep="\t", header=None)
#df_viruses = pd.read_csv("/Users/enriquem.muro/Desktop/viruses/spNames_phages__stat_description.protein.uniprot_reference_proteome.txt", sep="\t", header=None)

if 0:
    pd.options.display.max_columns=None
    print(df_viruses.head(3))
    sys.exit()
virus_names=df_viruses[0].tolist()
print(virus_names)
pd.options.display.max_columns=None
#print(df_prot_stat_desc.head(2))
df_prot_stat_desc.loc[df_prot_stat_desc['species'].isin(virus_names), 'Human_virus'] = True
df_prot_stat_desc['Human_virus'] = df_prot_stat_desc['Human_virus'].astype(object)#
df_prot_stat_desc=df_prot_stat_desc.sort_values(by='Human_virus',ascending=True)

if 0:
    # show data types of all the columns
    print(df_prot_stat_desc.loc[df_prot_stat_desc['species'].isin(virus_names)])
    # Get a Series object containing the data type objects of each column of Dataframe.
    # Index of series is column name.
    dataTypeSeries = df_prot_stat_desc.dtypes
    print('Data type of each column of Dataframe :')
    print(dataTypeSeries)

#print(df_prot_stat_desc.groupby("Human_virus").count())
#sys.exit()
#
# sort by mean
df_hv = df_prot_stat_desc[df_prot_stat_desc['Human_virus'] == True]
df_hv = df_hv[["species", "mean"]]
df_hv.sort_values("mean", inplace=True)
if 1:
    pd.options.display.max_columns=None
    print("-->", len(df_hv))
    l_mean_hv=df_hv["mean"].tolist()
    print(l_mean_hv)
    import matplotlib.pyplot as plt
    plt.hist(l_mean_hv,100)
    plt.show()
    sys.exit()
    print(df_hv.tail(50))
    #my_list = df2plot["var"].to_list()
    #print(my_list) # check that there are no nan value
    sys.exit()

#
df2plot = df_prot_stat_desc.copy()
df2plot = df2plot[["superregnum", "count", "mean", "var","Human_virus"]]
df2plot = df2plot[df2plot.superregnum=="viruses"]
df2plot.dropna(subset=["mean", "var"], inplace=True)

df2plot = df2plot[df2plot["mean"]>0]
df2plot = df2plot[df2plot["var"]>0]

if 0:
    pd.options.display.max_columns=None
    print(df2plot.head(9))
    #my_list = df2plot["var"].to_list()
    #print(my_list)  # check that there are no nan value
    sys.exit()

if 1:
    #calculate best fit line
    slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(df2plot['mean']), np.log10(df2plot['var']))
    #format text
    txt = 'v = {:4.4} * m^{:4.4};   R^2= {:2.4f}'.format(10**intercept, slope, r_value**2)
else:
    txt="Nothing has been calculated"
print(txt)
#sys.exit()

p = (ggplot(df2plot, aes("mean", "var", color="Human_virus")) + geom_point(size=0.75) +
     geom_smooth(method="lm", color="green", size=0.2, span=.8) +
     labs(title="Reference proteomes (protein length)")
     + scale_x_log10(breaks=[10 ** power for power in range(6)],
                     limits=[min(df2plot["mean"].to_list())/2, 2*max(df2plot["mean"].to_list())]
                     #limits=[400, 600]
                   ) +
     scale_y_log10(breaks = [10**power for power in range(13)],
                   limits = [min(df2plot["var"].to_list())/2, 2*max(df2plot["var"].to_list())])
    ) +  theme(legend_position=(0.75, 0.3), legend_key_size=5, legend_background=element_rect(fill='grey', alpha=0.01)) + \
    annotate('text', x=0.05*max(df2plot["mean"].to_list()), y=0.9*max(df2plot["var"].to_list()),
             label=txt, size=7, color="black")
if 1:
    print(p)