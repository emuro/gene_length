import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd
import numpy as np
#from pandas.api.types import CategoricalDtype
from plotnine import *
import sys
#
#%matplotlib inline
#from scipy.stats import chisquare, chi2_contingency
#from statistics import print_results_chi2_contingency, print_results_2var_correlation


'''
def plot_dfColumn_histogram(df, col="length",
                            plot_title="", x_axis_label="",
                            bool_plot=1, plot_file_withPath=""): #  plot: 0, 1
    """
        input:
            -df containing the data to be plotted (sorted by number)
                length: lengths

            -bool_plot=1,          #  plot: 0, 1
            -plot_file_withPath="" # in the case bool_plot==1
    """
    if 0:
        pd.set_option('display.max_rows', None)
        print(df.head(10))
        print(df.dtypes())
        #
        print(df["contig"].value_counts())
        exit()
    #print(df_subSet["gene_id"].to_string(index=False))
    #print(df.length.dtypes)
    #exit()


    #df["length"])
    p=(ggplot(df)         # defining what data to use
     + aes(x=col)    # defining what variable to use
     + labs(x=x_axis_label, y="number of genes", title=plot_title) # legend title
     + geom_histogram(colour="black", alpha=0.5, fill="red") # defining the type of plot to use
     + scale_x_log10(breaks=[1,2,3,4,5,10,100,1000,10000,100000,1000000])   #xaxis: log10
     #+ scale_x_continuous(trans = 'log2')                                    #xaxis: log2
       )
     #p = p + theme_xkcd()

    if bool_plot and plot_file_withPath:
        p.save(plot_file_withPath, verbose=False) # save the plot
    else:
        print(p)

    return
'''


def plot_dfColumn_histogram_Log10x(df, col="log10ofLength", width_of_bin=0.2,
                                   count=1000, mean=4, std=2, min=1, max=6, # a ojo de buen cubero
                                   plot_title="", x_axis_label="",
                                   bool_plot=1, plot_file_withPath=""): #  plot: 0, 1
    """
        input:
            -df containing the data to be plotted (sorted by number)
                length: lengths

            -bool_plot=1,          #  plot: 0, 1
            -plot_file_withPath="" # in the case bool_plot==1
    """
    if 0:
        pd.set_option('display.max_rows', None)
        print(df.head(10))
        ###print(df.dtypes())
        #
        print(df["contig"].value_counts())
        sys.exit()
    #print(df_subSet["gene_id"].to_string(index=False))
    #print(df.length.dtypes)
    #exit()

    print("->", count, mean, std, min, max, sep="\t")


    x = np.linspace(0, 2*mean, count)
    print(x)
    pdf = np.exp(-(x-mean)**2 / (2 * std**2))/ (std * np.sqrt(2 * np.pi))
    print("--->",sum(pdf/(2*mean)))


    df_for_plotting = df[[col]].copy() #double square is necessary
    df_for_plotting["x"]   = x
    df_for_plotting["pdf"] = 1000*pdf  #by number of elements
    if 0:
        pd.options.display.max_columns = None
        pd.options.display.max_rows = None
        print(df_for_plotting.head(19))

    #df["length"])
    p=(ggplot(df_for_plotting)         # defining what data to use
     + aes(x=col)    # defining what variable to use
     + labs(x=x_axis_label, y="number of genes", title=plot_title) # legend title
     + geom_histogram(colour="black", alpha=0.5, fill="red", binwidth=width_of_bin) # defining the type of plot to use
     + geom_point(aes(x="x", y="pdf"), size=0.01, color="lightgreen")
     + scale_x_continuous(breaks=np.arange(0, 7, 1.), minor_breaks=np.arange(0, 7, .25), limits=[0, 7])
    )

    if bool_plot and plot_file_withPath:
        p.save(plot_file_withPath, verbose=False) # save the plot
    else:
        print(p)
    return

def plot_withMatplotLib_dfColumn_histogram_Log10x(data,
                                   count=1000, division="vertebrates", annotation_name="title",
                                   text__is_a_logNorm="no log-norm statistics",
                                   entity="gene", # gene or protein
                                   x_min=1, x_max=7, # a ojo de buen cubero
                                   # plot_title="", x_axis_label="",
                                   bool_plot=1, plot_file_withPath=""): #  plot: 0, 1
    """
        input:
            -data: list
                containing the data to be plotted
            -count: int
                number of genes
            -division: str
                the ensembl division of the species
            -text__is_a_logNorm: str
                text to be embedded in the figure
            -annotation_name: str
                the species
            -x_min, x_max: int
                min and max for the x_axis
            -bool_plot=1,          #  plot: 0, 1
            -plot_file_withPath="" # in the case bool_plot==1
    """

    mean, std = stats.distributions.norm.fit(data)
    # Plot the histogram.
    plt.hist(data, bins=50, density=True, alpha=0.6, color='b', edgecolor='black', linewidth=1.2)
    title = division + ": " + annotation_name
    plt.figtext(.5, .9, title, fontsize=8, ha='center') # little size because some species have a very long title.
    plt.xlabel(entity + " length")
    plt.ylabel("probability density function (pdf)")
    plt.text(0.55, 0.6,
             text__is_a_logNorm,
             transform=plt.gca().transAxes)
    plt.text(0.01, 0.7,
             "count = %d\nmean = %.2f\nstd = %.2f" % (count, mean, std), transform=plt.gca().transAxes)
    x = np.linspace(x_min, x_max, 700)
    fitted_data =  stats.distributions.norm.pdf(x, mean, std) # integral = 1
    plt.plot(x, fitted_data, 'r-')
    plt.vlines(mean, 0, max(fitted_data), colors='red', linestyles='dotted')
    plt.savefig(plot_file_withPath,dpi=300,bbox_inches='tight')
    if bool_plot:
        plt.savefig(plot_file_withPath, dpi=300, bbox_inches='tight')
    else:
        plt.show()
    plt.clf()
    plt.close()
    return

'''
def plot_sorted_histogram_of_biotypes_old(df, plot_title="plot_sorted_histogram_of_biotype",
                                      bool_plot=1,  #  plot: 0, 1
                                      plot_file_withPath=""):
    """
        input:
            -df containing the data to be plotted (sorted by number)
                biotype: names of gene types
                number: number of genes for each gene type

            -bool_plot=1,          #  plot: 0, 1
            -plot_file_withPath="" # in the case bool_plot==1
    """

    df.biotype = df.biotype.astype('category') # convert "biotype" column into "category" type.
    df['biotype'] = df['biotype'].cat.reorder_categories(df.biotype.tolist()) #sort...
    ##print(df.dtypes)

    p=(ggplot(df)         # defining what data to use
     + aes(x="biotype", y="number", fill="biotype")    # defining what variable to use
     + labs(x="type of gene", y="number of genes", title=plot_title, fill="") # legend title
     + geom_bar(stat="identity", width=1, colour="black", alpha=0.75, show_legend=True) # defining the type of plot to use
     + geom_text(aes(label="number"), size=7, va='bottom', ha='center')
    )
    p = p + theme(legend_direction="horizontal", legend_position=(0.65, .8),
                  legend_key_size=10, legend_background=element_rect(fill = 'grey', alpha=0.01))

    # place legend at top and grey axis lines
    p = p + theme(axis_line_y=element_line(color="black"),
                  axis_line_x=element_line(color="black"),
                  axis_ticks_major_x=element_line(size=0.25),
                  axis_ticks_minor_x=element_line(size=0.25),
                  axis_text_x=element_text(color="black", alpha=0.01),  # make it transparent :-)
                  axis_text_y=element_text()  # (color="black", size=13)
                  )
    p = p + theme(axis_text_x=element_text(angle=0, margin={'t': 0, 'r': 0}, size=7))
    if bool_plot and plot_file_withPath:
        p.save(plot_file_withPath, verbose=False) # save the plot
    else:
        print(p)

    return


def plot_sorted_histogram_of_biotypes(df, plot_title="plot_sorted_histogram_of_biotype",
                                      bool_plot=1,  #  [0:display|1:save to file]
                                      plot_file_withPath=""):
    """
        input:
            -df containing the data to be plotted (sorted by number)
                biotype: names of gene types
                number: number of genes for each gene type
        input and output:
            -bool_plot=1,          # 0|1: display|save file
            -plot_file_withPath="" # in the case bool_plot==1
    """

    df.biotype = df.biotype.astype('category')  # convert "biotype" column into "category" type.
    df['biotype'] = df['biotype'].cat.reorder_categories(df.biotype.tolist())  # sort...
    if 0:
        print(df.dtypes)

    p=(ggplot(df)         # defining what data to use
     + aes(x="biotype", y="number", fill="biotype")  # defining what variable to use
     + labs(x="type of gene", y="number of genes", title=plot_title, fill="")  # legend title
     + geom_bar(stat="identity", width=1, colour="black", alpha=0.75, show_legend=True) # the type of plot: bar
     + geom_text(aes(label="number"), size=4, color="black", va='bottom', ha='center')
    )
    p = p + theme(legend_direction="horizontal", legend_position=(0.54, .8),
                  legend_key_size=5, legend_text=element_text(size=3.75),
                  legend_background=element_rect(fill = 'grey', alpha=0.01))
    # place legend at top and grey axis lines
    p = p + theme(axis_line_y=element_line(color="black"),
                  axis_line_x=element_line(color="black"),
                  axis_ticks_major_x=element_line(size=0.25),
                  axis_ticks_minor_x=element_line(size=0.25),
                  axis_text_x=element_text(color="black", alpha=0.01), # make it transparent :-)
                  axis_text_y=element_text() #(color="black", size=13)
                  )
    p = p + theme(axis_text_x=element_text(angle=0, margin={'t': 0, 'r': 0}, size=7))
    if bool_plot and plot_file_withPath:
        p.save(plot_file_withPath, verbose=False, dpi=300) # save the plot to a file; dpi≤1000
    else:
        print(p)

    return


'''
def plot_singleComparison2Benford__from_Classes(distrib, benford, how_many_digits,
                                    bar_or_pointLine="point-line", #"bar" or "point-line" representation
                                    count_or_prob="prob",          #"count" or "prob" representation
                                    bool_plot_line=0,              # in the case  bar_or_pointLine="point-line"
                                                                 #  0 plot without the line (only points)
                                                                 #  1 plot the line
                                    bool_plot=1,                    #  plot: 0, 1
                                    plot_file_withPath=""
                                  ):
    """
        input:
            receives 2 instances of DigitDistribution for plotting
            -distrib: the data observed
            -benford:  the benford law
            -how_many_digits: "1d" or "2d"
            -bar_or_pointLine: #"bar" or "point-line" representation
            -count_or_prob:    #"count" or "prob"
            -bool_plot_line    #when bar_or_pointLine="point-line"
                               #0 plot without the line (only points)
                               #  1 plot the line
            -bool_plot=1,      #  plot: 0, 1
            -plot_file_withPath="" # in the case bool_plot==1

        output
            The chi-sqrt variables:
                - x_sqrt
                - pvalue
                - dof
    """

    if not distrib.number_of_lengths:
        print("ERROR...plot_singleComparison2Benford",
            "distrib.number_of_lengths is not inizialitated" )
        exit()
    else:
        numGenesLegend = str(distrib.number_of_lengths) + " genes"
        metadata = distrib.name.split(":") #metadata[2] is the species
        species = "("+metadata[2]+")"
        the_figure_title = metadata[0].title()+" "+ metadata[1]+" ("+metadata[2]+")"

    if how_many_digits=="1d":
        x_axis_breaks = range(1,10,1)
        x_axis_label = "first digit"
        message_x_axis = 2.
        message_y_axis = 1.75
    elif how_many_digits=="2d":
        x_axis_breaks = range(10,100,5)
        x_axis_label = "first two digits"
        message_x_axis = 20
        message_y_axis = 1.25
    else:
        print("ERROR...plot_singleComparison2Benford",
            "distrib.number_of_lengths is not inizialitated" )
        exit()

    dfB = benford.get_df__d_count_prob(how_many_digits)
    dfB["count"] = dfB.prob * distrib.number_of_lengths
    dfB["count_int"] = dfB["count"].round()
    dfB["count_int"] = dfB["count_int"].astype(int)
    dfD = distrib.get_df__d_count_prob(how_many_digits)
    df = pd.concat([dfB, dfD])
    print(dfD)
    print(dfB)
    ##print(df.dtypes)

    #CHI-SQRT
    df_ctab = pd.DataFrame([ dfD["count"].to_list(),
                             dfB["count_int"].to_list()], columns=dfD["d"].tolist(), index=[dfD["biotype"][0], dfB["biotype"][0]])
    print_results_chi2_contingency(df_ctab, 0)
    x_sqrt, pvalue, dof, expected = chi2_contingency(df_ctab, correction=True)
    x_sqrt_tex = "{:1.2f}".format(x_sqrt)
    pvalue_tex = "{:.2e}".format(pvalue)
    text4plot_chisq = "X²:"+ x_sqrt_tex+ " df:"+ str(dof)+ "\np_value:" + pvalue_tex

    df_aux=pd.DataFrame(distrib.length, columns=["lengths"])
    ##print(df_aux.describe())

    if not bool_plot:
        return  x_sqrt, pvalue, dof
    if bar_or_pointLine=="point-line": # geom_point() + geom_line()
        if bool_plot_line==1:
            p = (ggplot(df, aes(x='d', y=count_or_prob, color="biotype", group="biotype"))
                + geom_point() +  geom_line()
                + labs(x=x_axis_label, y=count_or_prob, title=the_figure_title , color=numGenesLegend)
                + annotate(geom="text", x=message_x_axis , y=(df[count_or_prob].max()* message_y_axis/10), label=text4plot_chisq, color="black", size=10)
                + scale_x_continuous(breaks=x_axis_breaks)
             )
        elif bool_plot_line==0:
            p = (ggplot(df, aes(x='d', y=count_or_prob, color="biotype", group="biotype"))
                + geom_point()
                + labs(x=x_axis_label, y=count_or_prob, title=the_figure_title , color=numGenesLegend)
                + annotate(geom="text", x=message_x_axis , y=(df[count_or_prob].max()*message_y_axis/10), label=text4plot_chisq, color="black", size=10)
                + scale_x_continuous(breaks=x_axis_breaks)
             )
        else:
            print("ERROR...plot_singleComparison2Benford",
            "bool_plot_line is not boolean" )
            exit()
        p + labs(caption = "(based on data from ...)")
        p = p + theme(legend_direction="vertical", legend_position=(0.7, .8),
                      legend_title=element_text(size=9), legend_key_size=7,
                      legend_background=element_rect(fill='grey', alpha=0.01))
    elif bar_or_pointLine=="bar": # geom_bar()
        p = (ggplot(df, aes(x='d', y=count_or_prob, fill="biotype"))
            + geom_bar(position="dodge",stat="identity", colour="black", alpha=0.75)
            + labs(x=x_axis_label, y=count_or_prob, title=the_figure_title , fill=numGenesLegend)
            + scale_x_continuous(breaks=x_axis_breaks)
        )
        p = p + theme(legend_direction="vertical", legend_position=(0.7, .8),
                      legend_title=element_text(size=9), legend_key_size=7,
                      legend_background=element_rect(fill='grey', alpha=0.01))
    else:
        print("ERROR...plot_singleComparison2Benford", "bar_or_pointLine=",
            bar_or_pointLine)
        exit()
    ##print(p)
    if bool_plot and plot_file_withPath:
        p.save(plot_file_withPath, verbose=False) # save the plot
    return  x_sqrt, pvalue, dof

def plot_singleComparison2Benford(dfD, dfB, name_distribution_data, how_many_digits,
                                  bar_or_pointLine="point-line", #"bar" or "point-line" representation
                                  count_or_prob="prob",          #"count" or "prob" representation
                                  bool_plot_line=0,              # in the case  bar_or_pointLine="point-line"
                                                                 #  0 plot without the line (only points)
                                                                 #  1 plot the line
                                  bool_plot=1,                   #  plot: 0, 1
                                  plot_file_withPath=""
                                  ):
    """
        input:
            receives 2 instances of DigitDistribution for plotting
            -distrib: the data observed
            -benford:  the benford law
            -name_distribution_data: info about the distribution to compare with
            -how_many_digits: "1d" or "2d"
            -bar_or_pointLine: #"bar" or "point-line" representation
            -count_or_prob:    #"count" or "prob"
            -bool_plot_line    #when bar_or_pointLine="point-line"
                               #0 plot without the line (only points)
                               #  1 plot the line
            -bool_plot=1,      #  plot: 0, 1
            -plot_file_withPath="" # in the case bool_plot==1

        output
            The chi-sqrt variables:
                - x_sqrt
                - pvalue
                - dof
    """

    numGenesLegend = str(dfD["count"].sum()) + " genes"
    metadata = name_distribution_data.split(":") #metadata[2] is the species
    the_figure_title = metadata[0].title()+" "+ metadata[1]+" ("+metadata[2]+")"

    if how_many_digits=="1d":
        x_axis_breaks = range(1,10,1)
        x_axis_label = "first digit"
        message_x_axis = 2.
        message_y_axis = 1.75
    elif how_many_digits=="2d":
        x_axis_breaks = range(10,100,5)
        x_axis_label = "first two digits"
        message_x_axis = 20
        message_y_axis = 1.25
    else:
        print("ERROR...plot_singleComparison2Benford",
            "wrong number of digits"  )
        exit()

    dfB["count"] = dfB.prob * dfD["count"].sum()
    dfB["count_int"] = dfB["count"].round()
    dfB["count_int"] = dfB["count_int"].astype(int)
    df = pd.concat([dfB, dfD])
    if 0:
        print(dfD)
        print(dfB)

    #CHI-SQRT
    df_ctab = pd.DataFrame([ dfD["count"].to_list(),
                             dfB["count_int"].to_list()], columns=dfD["d"].tolist(), index=[dfD["biotype"][0], dfB["biotype"][0]])
    x_sqrt, pvalue, dof, expected = print_results_chi2_contingency(df_ctab, 0, 0) # 0, 1: expected and print
    x_sqrt_tex = "{:1.2f}".format(x_sqrt)
    pvalue_tex = "{:.2e}".format(pvalue)
    text4plot_chisq = "X²:"+ x_sqrt_tex+ " df:"+ str(dof)+ "\npvalue:" + pvalue_tex

    if not bool_plot:
        return  x_sqrt, pvalue, dof
    if bar_or_pointLine=="point-line": # geom_point() + geom_line()
        if bool_plot_line==1:
            p = (ggplot(df, aes(x='d', y=count_or_prob, color="biotype", group="biotype"))
                + geom_point() +  geom_line()
                + labs(x=x_axis_label, y=count_or_prob, title=the_figure_title , color=numGenesLegend)
                + annotate(geom="text", x=message_x_axis , y=(df[count_or_prob].max()* message_y_axis/10), label=text4plot_chisq, color="black", size=10)
                + scale_x_continuous(breaks=x_axis_breaks)
             )
        elif bool_plot_line==0:
            p = (ggplot(df, aes(x='d', y=count_or_prob, color="biotype", group="biotype"))
                + geom_point()
                + labs(x=x_axis_label, y=count_or_prob, title=the_figure_title , color=numGenesLegend)
                + annotate(geom="text", x=message_x_axis , y=(df[count_or_prob].max()*message_y_axis/10), label=text4plot_chisq, color="black", size=10)
                + scale_x_continuous(breaks=x_axis_breaks)
             )
        else:
            print("ERROR...plot_singleComparison2Benford",
            "bool_plot_line is not boolean" )
            exit()
        p + labs(caption = "(based on data from ...)")
        p = p + theme(legend_direction="vertical", legend_position=(0.7, .8),
                      legend_title=element_text(size=9), legend_key_size=7,
                      legend_background=element_rect(fill='grey', alpha=0.01))
    elif bar_or_pointLine=="bar": # geom_bar()
        p = (ggplot(df, aes(x='d', y=count_or_prob, fill="biotype"))
            + geom_bar(position="dodge",stat="identity", colour="black", alpha=0.75)
            + labs(x=x_axis_label, y=count_or_prob, title=the_figure_title , fill=numGenesLegend)
            + scale_x_continuous(breaks=x_axis_breaks)
        )
        p = p + theme(legend_direction="vertical", legend_position=(0.7, .8),
                      legend_title=element_text(size=9), legend_key_size=7,
                      legend_background=element_rect(fill='grey', alpha=0.01))
    else:
        print("ERROR...plot_singleComparison2Benford", "bar_or_pointLine=",
            bar_or_pointLine)
        exit()
    if bool_plot and plot_file_withPath:
        p.save(plot_file_withPath, verbose = False) # save the plot
    else:
        print(p)
    return  x_sqrt, pvalue, dof


def plot_2var_correlation(dfD, dfB, name_distribution_data, how_many_digits,
                          bool_plot=0,                   #  0 plot without the line (only points)
                          plot_file_withPath=""          #  1 plot the line
                        ):
    """
        input:
            receives 2 instances of DigitDistribution for plotting
            -dfD: data frame, the data observed
            -dfB: dataframe, the benford law
            -how_many_digits: "1d" or "2d"
            -bool_plot=1,      #  plot: 0, 1
            -plot_file_withPath="" # in the case bool_plot==1
        output
            The pearson variables:
                - pearson's correlation
                - pvalue
    """
    numGenesLegend = str(dfD["count"].sum()) + " genes"
    metadata = name_distribution_data.split(":") #metadata[0] is the db
                                                 #metadata[1] db version
                                                 #metadata[2] is the species
                                                 #metadata[3] is the genetype
    the_figure_title = metadata[0].title()+" "+ metadata[1]+" ("+metadata[2]+")"
    if how_many_digits=="1d":
        pass
    elif how_many_digits=="2d":
        pass
    else:
        print("plot_2var_correlation",
            "wrong number of digits" )
        exit()

    dfB["count"] = dfB.prob * dfD["count"].sum()
    dfB["count_int"] = dfB["count"].round()
    dfB["count_int"] = dfB["count_int"].astype(int)
    if 0:
        print(dfD)
        print(dfB)

    #Pearson correlation
    df_ctab = pd.DataFrame([ dfD["count"].to_list(),
                             dfB["count"].to_list()], columns=dfD["d"].tolist(), index=[dfD["biotype"][0], dfB["biotype"][0]])
    if 0:
        print(df_ctab.T)
    df_ctab_trans = df_ctab.T
    c = df_ctab_trans.columns
    v1 = df_ctab_trans[c[0]].to_list()
    v2 = df_ctab_trans[c[1]].to_list()
    pearson_corr, pvalue_corr = print_results_2var_correlation(v1, v2, 0) # 0, 1 to print in the screen
    text4plot_pearson = "Pearson's correlation: "+str("{:.6f}".format(pearson_corr))+"\npvalue: "+ str("{:.2e}".format(pvalue_corr))
    #
    plotdf = pd.DataFrame()
    plotdf["Benford"] = dfB["count"]
    plotdf["observation"] = dfD["count"]
    max_val = max(plotdf["Benford"].to_list() + plotdf["observation"].to_list())

    p = ( ggplot(plotdf, aes(y="Benford", x="observation")) # we could also plot the prob instead of freq (prob.d)
            + xlab("Num. of genes ("+str(plotdf["observation"].sum())+" "+metadata[3]+")")
            + ylab("Num. of genes (Benford)")
            #+ xlim(0, max_val) + ylim(0, max_val)
            + labs(caption="", title=the_figure_title)  # same format than other plots, other (caption = database, title=mainTitle)
            + geom_smooth(method='lm', span=0.3, size = 0.2, colour='lightblue', se=True, alpha=0.5 )
            + geom_abline(slope=1, colour='lightgreen', size = 1.0, linetype = "solid") # compare to benford #dashed
            + geom_point(size = 1.0)
            + annotate(geom="text", x=0.35*max_val, y=0.95*max_val, label=text4plot_pearson, color="black", size=10)
        )
    if bool_plot and plot_file_withPath:
        p.save(plot_file_withPath, verbose = False) # save the plot
    else:
      print(p)
    return pearson_corr, pvalue_corr


def prepare_plot_compare2Benford(distrib, benford, how_many_digits):
    """
        input:
            receives 2 instances of DigitDistribution for plotting
            -distrib: the data observed
            -benford:  the benford law
            -how_many_digits: "1d" or "2d"
        output
            two df :
                - dfD: data frame for the Data
                - dfB: data frame for the Benford
    """

    if not distrib.number_of_lengths:
        print("ERROR...plot_2var_correlation",
            "distrib.number_of_lengths is not inizialitated" )
        exit()

    if how_many_digits not in ["1d", "2d"]:
        print("ERROR...plot_singleComparison2Benford",
            "distrib.number_of_lengths is not inizialitated" )
        exit()

    dfB = benford.get_df__d_count_prob(how_many_digits)
    dfB["count"] = dfB.prob * distrib.number_of_lengths
    dfB["count_int"] = dfB["count"].round()
    dfB["count_int"] = dfB["count_int"].astype(int)
    dfD = distrib.get_df__d_count_prob(how_many_digits)
    if 0:
        print(dfD)
        print(dfB)
    return dfD, dfB
