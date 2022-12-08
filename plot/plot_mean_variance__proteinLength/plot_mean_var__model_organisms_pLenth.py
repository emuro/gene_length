from lib_analysis import constants_analysis as c
import pandas as pd
import sys
import numpy as np

#import plotnine as p9
from plotnine import *
print(ggplot.__doc__)
from scipy import stats


df = pd.DataFrame()
if 1:
    stat_file = c.LOG_SOME_STATISTICS_PROTEIN_FILE
else:
    stat_file = c.LOG_SOME_STATISTICS_REFERENCE_PROTEOME_FILE

if 0:
    print(stat_file)
df_prot_stat_desc = pd.read_csv(stat_file, sep="\t")
if 0:
    pd.options.display.max_columns=None
    print(df_prot_stat_desc.head(9))
    sys.exit()

df2plot = df_prot_stat_desc.copy()

#mean vs. var
#calculate best fit line
slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(df2plot['mean']), np.log10(df2plot['var']))
#format text
txt = 'v = {:4.4} * m^{:4.4};   R^2= {:2.4f}'.format(10**intercept, slope, r_value**2)
print(txt)
p = (ggplot(df2plot, aes("mean", "var", color="species")) + geom_point(size=0.75) +
     geom_smooth(method="lm", color="green", size=0.2, span=.8)
     + labs(title="Proteins (some model organisms)")
     + scale_x_log10(breaks=[10 ** power for power in range(6)],
                   limits=[min(df2plot["mean"].to_list())/2, 2*max(df2plot["mean"].to_list())]
                   ) +
     scale_y_log10(breaks = [10**power for power in range(13)],
                   limits = [min(df2plot["var"].to_list())/2, 2*max(df2plot["var"].to_list())])
    ) +  theme(legend_position=(0.75,0.3), legend_key_size=5, legend_background=element_rect(fill='grey', alpha=0.01)) + \
    annotate('text', x=0.45*max(df2plot["mean"].to_list()), y=0.95*max(df2plot["var"].to_list()),
             label=txt, size=7, color="black")
if 1:
    print(p)

#count vs. mean
#calculate best fit line
slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(df2plot['count']), np.log10(df2plot['mean']))
#format text
txt = 'v = {:4.4} * m^{:4.4};   R^2= {:2.4f}'.format(10**intercept, slope, r_value**2)
print(txt)
p2 = (ggplot(df2plot, aes("count", "mean", color="species")) + geom_point(size=0.75) +
     geom_smooth(method="lm", color="green", size=0.2, span=.8)
     + labs(title="Proteins (some model organisms)")
    ) +  theme(legend_position=(0.75,0.725), legend_key_size=5, legend_background=element_rect(fill='grey', alpha=0.01)) + \
    annotate('text', x=0.45*max(df2plot["count"].to_list()), y=1.5*max(df2plot["mean"].to_list()),
             label=txt, size=7, color="black")
if 0:
    print(p2) # for some reason, I have to print my plot

#count vs. mean
#calculate best fit line
slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(df2plot['count']), np.log10(df2plot['var']))
#format text
txt = 'v = {:4.4} * m^{:4.4};   R^2= {:2.4f}'.format(10**intercept, slope, r_value**2)
print(txt)
p3 = (ggplot(df2plot, aes("count", "var", color="species")) + geom_point(size=0.75) +
     geom_smooth(method="lm", color="green", size=0.2, span=.8)
     + labs(title="Proteins (some model organisms)")
    ) +  theme(legend_position=(0.3,0.7), legend_key_size=0.5, legend_background=element_rect(fill='grey', alpha=0.01)) + \
    annotate('text', x=0.6*max(df2plot["count"].to_list()), y=1.6*max(df2plot["var"].to_list()),
             label=txt, size=7, color="black")
if 0:
    print(p3) # for some reason, I have to print my plot