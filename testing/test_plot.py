from lib_analysis import constants_analysis as c
import pandas as pd
import numpy as np
import math

import sys
#import plotnine as p9
from plotnine import *
from scipy import stats



DIVISIONS = ["viruses", "bacteria", "protist", "fungi", "plants", "metazoa", "vertebrates"]
#DIVISIONS.reverse()
DBS = ["ensemblgenomes", "ensemblgenomes", "ensemblgenomes", "ensemblgenomes", "ensemblgenomes", "ensemblgenomes",
       "ensembl"] # for viruses, bateria, ...respectively
#DBS.reverse()
ENSEMBL_VERSIONS = ["101", "49", "49", "49", "49", "49", "98"] # for viruses, bateria, ...respectively
#ENSEMBL_VERSIONS.reverse()

df = pd.DataFrame()
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

df2plot = df.copy()

#calculate best fit line
slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(df2plot['mean']), np.log10(df2plot['var']))
#format text
txt = 'v = {:4.4} * m^{:4.4};   R^2= {:2.4f}'.format(10**intercept, slope, r_value**2)
p = (ggplot(df2plot, aes("mean", "var", color="Division")) + geom_point(size=0.75) +
     geom_smooth(method="lm", color="green", size=0.5, span=.8)
     + scale_x_log10(breaks=[10 ** power for power in range(6)],
                   limits=[min(df2plot["mean"].to_list())/2, 2*max(df2plot["mean"].to_list())]
                   ) +
     scale_y_log10(breaks = [10**power for power in range(13)],
                   limits = [min(df2plot["var"].to_list())/2, 2*max(df2plot["var"].to_list())])
    ) +  theme(legend_position=(0.8,0.3), legend_key_size=5, legend_background=element_rect(fill='grey', alpha=0.01)) + \
    annotate('text', x=0.25*max(df2plot["mean"].to_list()), y=0.3*max(df2plot["var"].to_list()),
             label=txt, size=7, color="darkgreen")

print(p) # for some reason, I have to print my plot