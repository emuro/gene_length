from lib import constants as c,use_pyensembl as pyens
from old_benford import benford_mod
#
import pandas as pd
import numpy as np
from time import time

# output data path directory
out_data_path=c.OUTPUT_PROJECT_PATH+c.DATA_PATH_NAME+"/"

#
# get all the info from all the species: info obtained from the indexing
df_species = pd.read_csv(c.LOG_ANNOTATION_INDEXED_SPECIES_FILE, sep="\t")
# eliminate undesired species
df_species = df_species[~df_species.species.isin(c.ANNOTATION_NAMES_TO_EXCLUDE)]
if 0:
    pd.options.display.max_columns = None
    print(df_species.head(2))
if 0:
    print(df_species.columns)
df_species.loc[:, ["count", "mean", "std", "min", "25perc", "50perc", "75perc", "max",
                   "log10_mean", "log10_std", "log10_min", "log10_25perc", "log10_50perc", "log10_75perc", "log10_max",
                   "log_mean", "log_std", "log_min", "log_25perc", "log_50perc", "log_75perc", "log_max"]] = "NaN"


dict_of_all_df_species = {}
list_of_all_biotypes = []
# iterate through each row (species)
df_biotypes_per_specie = pd.DataFrame()
#
# iterate through each row (species)
for i in range(len(df_species)):
    start_time_species=time()
    annotation_name = df_species.loc[i, "species"]
    reference_name  = df_species.loc[i, "assembly"]
    annotation_path = df_species.loc[i, "path"]
    annotation_file = df_species.loc[i, "file"]
    if c.BOOL_CHECK_SOME_SPECIES: #for the sake of a unique species analysis
        if annotation_name not in c.ANNOTATION_NAMES_TO_CHECK:
            continue
    print(df_species.loc[i, "species"], df_species.loc[i, "assembly"],
          df_species.loc[i,"ensembl_version"], df_species.loc[i, "path"],
          df_species.loc[i, "file"])

    # get ensembl object of gene annotations from the pyemsembl file that
    # I have previously indexed
    data = pyens.index_db_of_geneAnnotation_from_gtfFile(reference_name, annotation_name,
                                                         annotation_path+annotation_file)
    # get the species
    genes = pyens.retrieve_geneAnnotationObject_fromIndexedData(data)
    print("There are", len(genes), "species in total (species: raw data)")
    if 0:
        print(genes[2])

    # get df annotation of all species (calculate and include the length)
    # Also count the species by biotype (gene type)
    df_all_genes = pyens.geneObjectList_to_df_sorted_by_biotype_and_geneLength(genes)
    if 0:
        print(df_all_genes.head(10))
    df_sorted_count_biotype=pyens.count_genes_by_biotype_givenDf(df_all_genes,
                                                                 # LUEGO, despues de filtrar species!!!!!!!!!...
                                                                 c.MIN_NUM_OF_GENES_PER_BIOTYPE)
    sortedByCounts_biotype=df_sorted_count_biotype["biotype"].to_list()
    if 0:
        print(df_sorted_count_biotype)
        print(sortedByCounts_biotype)
    # Note this column has 4 rows of values:
    df_add_biotype_per_specie__aux = pd.DataFrame({
        annotation_name: sortedByCounts_biotype
    })
    df_biotypes_per_specie=pd.concat([df_biotypes_per_specie,df_add_biotype_per_specie__aux],ignore_index=False,
                                     axis=1)
    if 0:
        print(df_biotypes_per_specie)
    #
    # filter by gene types with more than c.MIN_NUM_OF_GENES_PER_BIOTYPE species
    df_all_genes = df_all_genes[(df_all_genes["biotype"].isin(sortedByCounts_biotype))]
    print("\t...after filtering by MIN_NUM_OF_GENES_PER_BIOTYPE (",
          c.MIN_NUM_OF_GENES_PER_BIOTYPE, "): there are ", len(df_all_genes),
          " species in total (df_all_genes)", sep="")

    #
    # iterate throughout each gene type
    for gene_type in sortedByCounts_biotype:
        if 1 and (gene_type not in c.CHECK_GENE_TYPE):  # just for several gene types
            continue

        if (gene_type not in list_of_all_biotypes):
            list_of_all_biotypes.append(gene_type)
            dict_of_all_df_species[gene_type]=df_species.copy()

        df_genes=df_all_genes[(df_all_genes["biotype"]==gene_type)]
        # the next line should eliminate whitespaces at the end
        # taken from https://stackoverflow.com/questions/33788913/pythonic-efficient-way-to-strip-whitespace-from-every-pandas-data-frame-cell-tha/33789292
        #df_genes=df_genes.apply(lambda x: x.str.strip() if x.dtype=="object" else x)

        print("There are ",len(df_genes)," species of this type (",gene_type,")",sep="")


        #
        # Calculate the difference of length to the closer gene
        #######################################################
        if 1:
            df_genes=df_genes.reset_index(drop=True) # I do not really need to reset the index
            #
            # get d list: the diff in length to any other gene
            l1=df_genes['length'].to_list()
            l=np.array(l1)
            lprev=np.roll(l,1)
            lnext=np.roll(l,-1)
            d=lnext-l  # distance with the next
            d[-1]=999999999
            dprev=l-lprev  # distance with the prev
            dprev[0]=999999999
            d[np.less(dprev,d)]=dprev[np.less(dprev,d)] #d is <class 'numpy.ndarray'>

            #The nex inserts avoid the SettingWithCopyWarning warning
            df_genes.insert(len(df_genes.columns),"diffLength",d)
            df_genes.insert(len(df_genes.columns),"log10ofLength",np.log10(df_genes.length))
            df_genes.insert(len(df_genes.columns),"logofLength",np.log(df_genes.length))
            if 0:
                pd.set_option('display.max_columns',None)
                print(df_genes.head(5))
                print(df_genes.tail(5))

            start_time=time()
            out_file=c.DB+"_"+c.ENSEMBL_VERSION+"."+annotation_name+"."+gene_type
            benford_mod.save_columns_givenDf(df_genes,
                                             ["biotype","length","diffLength","log10ofLength","logofLength","contig","start","end",
                                              "gene_id","gene_name"],
                                             1,out_data_path,out_file)
            elapsed_time=time()-start_time
            print(annotation_name," (", gene_type,
                  ")\n\tsaving length and contig files:\tElapsed time: %.10f seconds."%elapsed_time, sep="")
            del df_genes["diffLength"]  # I do not know if this is necessary ...but just in case
            del d

            #
            # Stat description
            ##################
            stat_describe_length=df_genes[["length","log10ofLength","logofLength"]].describe()
            if 0:
                print(stat_describe_length)
                print(stat_describe_length["length"]["std"])
            dict_of_all_df_species[gene_type].loc[i,"count"]=int(stat_describe_length["length"]["count"])
            #
            dict_of_all_df_species[gene_type].loc[i,"mean"]  =stat_describe_length["length"]["mean"]
            dict_of_all_df_species[gene_type].loc[i,"std"]   =stat_describe_length["length"]["std"]
            dict_of_all_df_species[gene_type].loc[i,"min"]   =int(stat_describe_length["length"]["min"])
            dict_of_all_df_species[gene_type].loc[i,"25perc"]=stat_describe_length["length"]["25%"]
            dict_of_all_df_species[gene_type].loc[i,"50perc"]=stat_describe_length["length"]["50%"]
            dict_of_all_df_species[gene_type].loc[i,"75perc"]=stat_describe_length["length"]["75%"]
            dict_of_all_df_species[gene_type].loc[i,"max"]   =int(stat_describe_length["length"]["max"])
            #
            dict_of_all_df_species[gene_type].loc[i,"log10_mean"]  =stat_describe_length["log10ofLength"]["mean"]
            dict_of_all_df_species[gene_type].loc[i,"log10_std"]   =stat_describe_length["log10ofLength"]["std"]
            dict_of_all_df_species[gene_type].loc[i,"log10_min"]   =stat_describe_length["log10ofLength"]["min"]
            dict_of_all_df_species[gene_type].loc[i,"log10_25perc"]=stat_describe_length["log10ofLength"]["25%"]
            dict_of_all_df_species[gene_type].loc[i,"log10_50perc"]=stat_describe_length["log10ofLength"]["50%"]
            dict_of_all_df_species[gene_type].loc[i,"log10_75perc"]=stat_describe_length["log10ofLength"]["75%"]
            dict_of_all_df_species[gene_type].loc[i,"log10_max"]   =stat_describe_length["log10ofLength"]["max"]
            #
            dict_of_all_df_species[gene_type].loc[i,"log_mean"]  =stat_describe_length["logofLength"]["mean"]
            dict_of_all_df_species[gene_type].loc[i,"log_std"]   =stat_describe_length["logofLength"]["std"]
            dict_of_all_df_species[gene_type].loc[i,"log_min"]   =stat_describe_length["logofLength"]["min"]
            dict_of_all_df_species[gene_type].loc[i,"log_25perc"]=stat_describe_length["logofLength"]["25%"]
            dict_of_all_df_species[gene_type].loc[i,"log_50perc"]=stat_describe_length["logofLength"]["50%"]
            dict_of_all_df_species[gene_type].loc[i,"log_75perc"]=stat_describe_length["logofLength"]["75%"]
            dict_of_all_df_species[gene_type].loc[i,"log_max"]   =stat_describe_length["logofLength"]["max"]

    elapsed_time=time()-start_time_species
    print(annotation_name," (total time): %.10f seconds.\n"%elapsed_time, sep="")

for b in list_of_all_biotypes:
    output_file = c.DB+"_"+c.ENSEMBL_VERSION+".species."+c.DATA_PATH_NAME+"."+b+".tsv"
    dict_of_all_df_species[b].to_csv(out_data_path+output_file, index=False, sep="\t") #log, only the first indexing
#
geneTypePerSpecie_out_file = c.DB+"_"+c.ENSEMBL_VERSION+"."+c.GENE_TYPE_PER_SPECIE+".tsv"
df_biotypes_per_specie.to_csv(out_data_path+geneTypePerSpecie_out_file,
                              index=False, sep="\t", na_rep='None')
if 0:
    print(df_biotypes_per_specie)