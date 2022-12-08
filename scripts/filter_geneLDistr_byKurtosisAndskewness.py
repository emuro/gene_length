from lib_analysis import constants_analysis as ca
import pandas as pd
import numpy as np
import sys
from shutil import copyfile



# i.e /Volumes/Wes/results/geneLength/outputInputFiles/analysis/some_statistics/stat_description/all.logNorm_stat.description.ensembl.tsv
all_logNorm_out_file = ca.OUTPUT_INPUT_FILES_PATH + ca.SOME_STATISTICS_PATH_NAME + \
                       "/all.logNormStat.description.ensembl.tsv"
df_stat = pd.read_csv(all_logNorm_out_file,sep="\t")
df_stat["new_kurtosis"] = df_stat["kurtosis"].abs()
df_stat["new_skew"] = abs(df_stat["skew"])
df_stat.sort_values(by=["new_kurtosis","new_skew"],ascending=True,inplace=True)
#
# condition to filter based in how good fits the log-norm
if 0:
    cond =df_stat["new_kurtosis"]>=0.025
    df_stat.drop(df_stat[cond].index,inplace=True)
    cond2 =df_stat["new_skew"]>=0.025
    df_stat.drop(df_stat[cond2].index,inplace=True)
    print(df_stat.shape)
#
if 0:
    print(all_logNorm_out_file)
    pd.options.display.max_columns=None
    print(df_stat.shape)
    print(df_stat.head(13))
    print(df_stat.columns)
    print(str(len(pd.unique(df_stat['species'])))+" diff species")
    sys.exit()

#d = dict(tuple(df_stat.groupby(df_stat['new_kurtosis'].diff().gt(0.05).cumsum())))
#pd.options.display.max_columns=None
#print (d[0][5])


bins=np.arange(0,1.0,0.05)
ind=np.digitize(df_stat['new_kurtosis'],bins)
# df_stat.groupby(ind).head()
#
g = df_stat.groupby(ind).groups
#print(df_stat.groupby(ind).groups)
for i in g.keys():
    print (i)
    print(g[i])
    for j in g[i]:
        png_in_file = df_stat.loc[j,"genes_file"].replace(".genes.",".geneLength_distrib.")
        png_in_file = png_in_file.replace(".tsv",".png")
        png_in_file_with_path = ca.OUT_DATA_LOCAL_PATH_ROOT + df_stat.loc[j,"trunk_genes_path"] + png_in_file

        png_out_file = "S"+str(i)+"_"+str(df_stat.loc[j,"kurtosis"])+"_"+str(df_stat.loc[j,"skew"])+"_"+str(df_stat.loc[j,"division"])+"_"
        png_out_file += df_stat.loc[j,"genes_file"].replace(".genes.",".geneLength_distrib.")
        png_out_file = png_out_file.replace(".tsv",".png")
        png_out_file_with_path = "/Users/enriquem.muro/tmp/gene_distr/" + png_out_file

        #print(str(i) + "\t" + str(df_stat.loc[j, "new_kurtosis"]) + "\t" + png_in_file_with_path + "\t" + str(df_stat.loc[j, "division"]))
        #print(str(i)+"\t"+str(df_stat.loc[j,"trunk_genes_path"])+"\t"+str(df_stat.loc[j,"kurtosis"])+"\t"+str(df_stat.loc[j,"skew"]))
        if 1:
            print(png_in_file_with_path)
            print(png_out_file_with_path)
        copyfile(png_in_file_with_path, png_out_file_with_path)