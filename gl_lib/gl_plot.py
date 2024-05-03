# python3
# ################################################################## #
# gl_plot.py (C) 2023 Mainz.
# Author: Enrique M. Muro
# ################################################################## #
#
# --------------------------------------------------------------------
# Project: gene_length
#
# ################################################################## #

import numpy as np
import pandas as pd
from scipy import stats
from plotnine import *
if 1:
    import warnings
    warnings.filterwarnings('ignore')

import sys
sys.path.append('./gl_lib')
import gl_constants as c

def plot_taylor_genes (df2plot, col_x, col_y, x_lab, y_lab, title, legends_by,\
     bool_show_regression): 
    #Calculate best fit line
    slope, intercept, r_value, p_value, std_err = \
        stats.linregress(np.log10(df2plot[col_x]),np.log10(df2plot[col_y]))
    #Format the regression text
    if 0: 
        print("v = {:4.4} * m^{:4.4};   R^2= {:2.4f}".format(10**intercept, slope, r_value**2))
    if bool_show_regression:
        txt = '$\sigma^{2} = ' + '{:4.2} '.format(10**intercept)  + ' \t ' +\
             '\mu^{' + '{:4.3}'.format(slope) + '}' + ';\tR^{2} = ' + '{:2.2f}$'.format(r_value**2)
    else:
        txt = ''
   
    p = (   
        ggplot(df2plot, aes(col_x, col_y, color=legends_by)) + 
        geom_point(size=0.15) + #, alpha=0.4 + # color, fill )
        geom_smooth(method="lm", color="green", size=0.25, span=.8) +
        labs(title=title, x=x_lab, y=y_lab) 
        + scale_color_manual(values=c.COLOR_ORG_GROUPS) 
        + labs(color='Group of organisms') # legend title
        + scale_x_log10(breaks=[10 ** power for power in range(6)],
          limits=[min(df2plot[col_x].to_list())/2, 2*max(df2plot[col_x].to_list())]) 
        + scale_y_log10(breaks = [10**power for power in range(13)], 
          limits = [min(df2plot[col_y].to_list())/2,2*max(df2plot[col_y].to_list())])#, labels=scientific_format(digits=2)
    ) + theme(legend_position=(0.95,0.1), legend_key_size=2, \
        legend_background=element_rect(fill='grey', alpha=0.01)) + \
            annotate('text', x=0.01*max(df2plot[col_x].to_list()), \
                y=0.65*max(df2plot[col_y].to_list()), label=txt,size=9,color="black")
    p.show()


def plot_taylor_proteins(df2plot, col_x, col_y, x_lab, y_lab, title, \
    legends_by, bool_show_regression): 

    # Calculate best fit line
    slope, intercept, r_value, p_value, std_err = \
        stats.linregress(np.log10(df2plot[col_x]),np.log10(df2plot[col_y]))
    # Format the regression text
    print('v = {:4.4} * m^{:4.4};   R^2= {:2.4f}'.format(10**intercept, slope, r_value**2))
    if bool_show_regression:
        txt = '$\sigma^{2} = ' + '{:4.2} '.format(10**intercept)  + '\t \mu^{' + \
            '{:4.3}'.format(slope) + '}' + ';\tR^{2} = ' + '{:2.2f}$'.format(r_value**2)
    else:
        txt = ''
        
    # limits and tick-breaks
    x_limits=[min(df2plot[col_x].to_list())/1.3, 1.3*max(df2plot[col_x].to_list())]
    y_limits=[min(df2plot[col_y].to_list())/1.3, 1.3*max(df2plot[col_y].to_list())]
    x_text_factor = 0.05
    y_text_factor = 0.95

    palette_colors = c.COLOR_ORG_GROUPS
    if legends_by == "Three-domain system":
        palette_colors = c.COLOR_KINGDOMS
    p = (   
        ggplot(df2plot, aes(col_x, col_y, color=legends_by)) + geom_point(size=0.1)
        + geom_smooth(method="lm", color="green", size=0.25, span=.8)
        + labs(title=title, x=x_lab, y=y_lab) 
        + scale_color_manual(values=palette_colors) # + scale_color_brewer() '#5C2D91'
        + labs(color=legends_by) # legend title
        + scale_x_log10(breaks=[10 ** power for power in range(2, 4)], limits=x_limits) #, labels=scientific_format(digits=2)
        + scale_y_log10(breaks = [10**power for power in range(3, 7)], limits=y_limits)
    ) + theme(legend_position=(0.95,0.1), legend_key_size=5, 
              legend_background=element_rect(fill='grey', alpha=0.01)) + \
                annotate('text', x=x_text_factor*max(df2plot[col_x].to_list()), \
                    y=y_text_factor*max(df2plot[col_y].to_list()), label=txt,size=9,color="black")
    p.show()


def plot_momentum (df2plot, col_x, col_y, x_lab, y_lab, title, \
    legends_by, legends_label, bool_show_regression, bool_prots=False): 
    #Calculate best fit line
    slope, intercept, r_value, p_value, std_err = \
        stats.linregress(np.log10(df2plot[col_x]),np.log10(df2plot[col_y]))
    #Format the regression text
    if 0: 
        print("v = {:4.4} * m^{:4.4} - m^{2};   \
            R^2= {:2.4f}".format(10**intercept, slope, r_value**2))
    if bool_show_regression:
        txt = '$\sigma^{2} = ' + '{:4.2} '.format(10**intercept)  + ' \t ' + \
            '\mu^{' + '{:4.3}'.format(slope) + '} - \mu^{2}' + ';\tR^{2} = ' + \
                '{:2.2f}$'.format(r_value**2)
    else:
        txt = ''
   

    # some plots parameter for genes or prots
    x_limits=[min(df2plot[col_x].to_list())/2, 2*max(df2plot[col_x].to_list())]
    y_limits=[min(df2plot[col_y].to_list())/2, 2*max(df2plot[col_y].to_list())]
    x_text_factor = 0.01
    y_text_factor = 0.65
    if bool_prots:
        # limits and tick-breaks
        x_limits=[min(df2plot[col_x].to_list())/1.3, 1.3*max(df2plot[col_x].to_list())]
        y_limits=[min(df2plot[col_y].to_list())/1.3, 1.3*max(df2plot[col_y].to_list())]
        x_text_factor = 0.05
        y_text_factor = 0.95
   
    palette_colors = c.COLOR_ORG_GROUPS
    if legends_by == "Three-domain system":
        palette_colors = c.COLOR_KINGDOMS
    p = (   
        ggplot(df2plot, aes(col_x, col_y, color=legends_by)) + 
        geom_point(size=0.15) + #, alpha=0.4) + # color, fill
        geom_smooth(method="lm", color="green", size=0.25, span=.8) +
        labs(title=title, x=x_lab, y=y_lab) 
                + scale_color_manual(values=palette_colors) 
        + labs(color=legends_label)
        + scale_x_log10(breaks=[10 ** power for power in range(6)],
          limits=x_limits) 
        + scale_y_log10(breaks = [10**power for power in range(13)], 
          limits =y_limits)#, labels=scientific_format(digits=2)
    ) + theme(legend_position=(0.95,0.1), legend_key_size=2, \
        legend_background=element_rect(fill='grey', alpha=0.01)) + \
        annotate('text', x=x_text_factor*max(df2plot[col_x].to_list()), \
            y=y_text_factor*max(df2plot[col_y].to_list()), label=txt,size=9,color="black")
    p.show()