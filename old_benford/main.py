from lib import rest_ensembl,constants as c

#import use_pyensembl as pyens
#import first_digit_distributions
#from plot import plot_dfColumn_histogram, plot_sorted_histogram_of_biotypes, prepare_plot_compare2Benford, plot_singleComparison2Benford, plot_2var_correlation
#import benford_mod
#import pandas as pd
#import numpy as np
#from statistics import print_results_Nigrini2000

DB = "ensembl"
ENSEMBL_VERSION = "98"

BOOL_CHECK_ONE_SPECIES   = c.BOOL_CHECK_ONE_SPECIES   #1 or 0
ANNOTATION_NAME_TO_CHECK = c.ANNOTATION_NAME_TO_CHECK #"caenorhabditis_elegans" #"danio_rerio" #"saccharomyces_cerevisiae" "drosophila_melanogaster" #"pan_paniscus" "mola_mola" #"mus_musculus" #"equus_caballus" #"anabas_testudineus" #"homo_sapiens" # species
REFERENCE_NAME_TO_CHECK  = c.REFERENCE_NAME_TO_CHECK  #"WBcel235" #"GRCz11" #"R64-1-1" #"panpan1.1" "ASM169857v1" #"GRCm38" #"EquCab3.0" #"fAnaTes1.1" #"GRCh38"
CHECK_GENE_TYPE          = c.CHECK_GENE_TYPE          #"[transcribed_processed_pseudogene", "protein_coding", "lncRNA"]

DB_LOCAL_PATH          = c.DB_LOCAL_PATH
PYCHARM_PROJECT_PATH   = c.PYCHARM_PROJECT_PATH

OUTPUT_PROJECT_PATH    = c.OUTPUT_PROJECT_PATH
#COMPARE_PATH_NAME   = "compare2Benford"
#CHI_SQRT_PATH_NAME  = "chiSqrt"
#CORR_PATH_NAME  = "correlation"
#GENE_TYPE_PATH_NAME = "geneTypeHist"
#GENE_TYPE_PER_SPECIE = "geneTypesPerSpecie"
#DATA_PATH_NAME      ="../../Desktop/benfordLaw/data"

MIN_NUM_OF_GENES_PER_BIOTYPE = c.MIN_NUM_OF_GENES_PER_BIOTYPE
LOG_ANNOTATION_INDEXED_SPECIES_FILE = c.LOG_ANNOTATION_INDEXED_SPECIES_FILE



# REST (diff. species)
######################
BOOL_REST_ENSEMBL_SPECIES = 1
if BOOL_REST_ENSEMBL_SPECIES:
    ensembl_release, number_of_species, df_species = rest_ensembl.get_REST_ensembl_species()
    print("\nREST ENSEMBL (rest_ensembl.get_REST_ensembl_species):")
    print("\tThe data corresponds to the release:", ensembl_release, "from ensembl")
    print("\tThere are", number_of_species, "different species in this release")
    print("\tdf_species sorted by taxonId:")
    rest_ensembl_species_file = PYCHARM_PROJECT_PATH + "rest_ensembl_" + str(ensembl_release) + "_species.tsv"
    df_species.to_csv(rest_ensembl_species_file, index = False, sep="\t")
    print("\t", df_species)
    exit()
dsagsdfas

#
# DB INDEXING
# NOTE: NEEDS TO BE RUN ONCE PER SESSION
# pyemsembl: Index the gene annotation for all the species given an
#            ensembl version. Once that it does it (~30 sec per species)
#            in the next session is extremely fast (secs for all species)
if 0: # index the encode annotation files for all the species
    BOOL_FIRST_INDEXING = 0
    df_species=use_pyensembl.indexData_for_all_species(ENSEMBL_VERSION, DB_LOCAL_PATH)
    if BOOL_FIRST_INDEXING:
        df_species.to_csv(LOG_ANNOTATION_INDEXED_SPECIES_FILE, index = False, sep="\t") #log, only the first indexing
    exit()


#
# get all the info from all the species: info obtained from the indexing
df_species = pd.read_csv(LOG_ANNOTATION_INDEXED_SPECIES_FILE, sep="\t")
if 0:
    pd.options.display.max_columns = None
    print(df_species.head(100))
if 0:
    print(df_species.columns)
df_species.loc[:, ["count", "mean", "std", "min", "25perc", "50perc", "75perc", "max",
                   "log10_mean", "log10_std", "log10_min", "log10_25perc", "log10_50perc", "log10_75perc", "log10_max",
                   "1d_xs_xsqrt", # new columns initialization
                   "1d_xs_pvalue", "1d_xs_dof", "2d_xs_xsqrt", "2d_xs_pvalue", "2d_xs_dof",
                   "1d_corr_pearson", "1d_corr_pvalue", "2d_corr_pearson", "2d_corr_pvalue",
                   "1d_nigrini2000Mad", "1d_nigrini2000_MAXad"]] = "NaN"

dict_of_all_df_species = {}
list_of_all_biotypes = []
#
# iterate through each row (species)
df_biotypes_per_specie = pd.DataFrame()
#
# iterate through each row (species)
for i in range(len(df_species)):
    annotation_name = df_species.loc[i, "species"]
    reference_name =  df_species.loc[i, "assembly"]
    annotation_file = df_species.loc[i, "path"]
    print(c.ANNOTATION_NAME_TO_CHECK.title())
    if BOOL_CHECK_ONE_SPECIES: #for the sake of a unique species analysis
        if annotation_name != c.ANNOTATION_NAME_TO_CHECK.title():
            continue
    print(df_species.loc[i, "species"], df_species.loc[i, "assembly"],
          df_species.loc[i, "ensembl_version"], df_species.loc[i, "path"])
    #
    # get ensembl object of gene annotations from the pyemsembl file that
    # I have previously indexed
    annotation_file = DB_LOCAL_PATH + "/pub/release-" + ENSEMBL_VERSION + "/gtf/" + annotation_name + "/" + annotation_file
    data = pyens.index_db_of_geneAnnotation_from_gtfFile(reference_name, annotation_name, annotation_file)
    # get the species
    genes = pyens.retrieve_geneAnnotationObject_fromIndexedData(data)
    print("There are", len(genes), "species in total (species: raw data)")
    if 0:
        print(genes[2])

    #
    # get df annotation of all species (calculate and include the length)
    # Also count the species by biotype (gene type)
    df_all_genes = pyens.geneObjectList_to_df_sorted_by_biotype_and_geneLength(genes)
    if 0:
        print(df_all_genes.head(10))
    df_sorted_count_biotype = pyens.count_genes_by_biotype_givenDf(df_all_genes, # LUEGO, despues de filtrar species!!!!!!!!!
                                                    MIN_NUM_OF_GENES_PER_BIOTYPE)
    sortedByCounts_biotype = df_sorted_count_biotype["biotype"].to_list()
    if 0:
        print(df_sorted_count_biotype)
    # Note this column has 4 rows of values:
    df_add_biotype_per_specie__aux = pd.DataFrame({
           annotation_name: sortedByCounts_biotype
    })
    df_biotypes_per_specie = pd.concat([df_biotypes_per_specie, df_add_biotype_per_specie__aux], ignore_index=False, axis=1)
    if 0:
        print (df_biotypes_per_specie)


    #
    # filter by gene types with more than MIN_NUM_OF_GENES_PER_BIOTYPE species
    df_all_genes = df_all_genes[(df_all_genes["biotype"].isin(sortedByCounts_biotype))]
    print("\t...after filtering by MIN_NUM_OF_GENES_PER_BIOTYPE (",
          MIN_NUM_OF_GENES_PER_BIOTYPE, "): there are ", len(df_all_genes),
          " species in total (df_all_genes)", sep="")


    # AQU'i hay que filtrarlos por otras cosas si fuera necesario  !!!!
    if 0: # for the sake of the chrms
        chroms, bad_chroms = pyens.get_standard_contigs_for_an_ensembl_version(data)
        ##print("chroms:", chroms)
        mask = df_all_genes.contig.isin(chroms)
        notMask = ~df_all_genes.contig.isin(chroms)
        df_genes = df_all_gene_onlyChrms = df_genes[mask]
        ##print(df_all_genes)

    # VER LOS DEL RANGO (parece que son olfatory receptors)
    if 0:
        df_all_genes_range = df_all_genes[(df_all_genes["length"] >= 915) & (df_all_genes["length"] <= 965)]
        print(df_all_genes_range)
        print(df_all_genes_range["gene_name"].to_string(index=False))
        exit()
    # ELIMINAR LOS DEL RANGO (parece que son olfatory receptors)
    if 0:
        df_all_genes  = df_all_genes[(df_all_genes["length"] <= 915) | (df_all_genes["length"] >= 965)]
        print(df_all_genes)
        #print(df_all_genes["gene_id"].to_string(index=False))
        #exit()
    # CAPTURAR LOS DEL RANGO (parece que son olfatory receptors)
    if 0:
        df_all_genes  = df_all_genes[(df_all_genes["length"] >= 915) & (df_all_genes["length"] <=965)]
        #print(df_all_genes)
        #exit()


    #
    # plot the gene type distribution
    BOOL_PLOT_GENETYPE_HIST = 1
    if BOOL_PLOT_GENETYPE_HIST:
        plot_path = OUTPUT_PROJECT_PATH+GENE_TYPE_PATH_NAME+"/"
        plot_file = DB+"_"+ENSEMBL_VERSION+"."+annotation_name+"."+GENE_TYPE_PATH_NAME+".png"
        plot_title = DB.capitalize() +" " + ENSEMBL_VERSION +" (" + annotation_name + ")"
        plot_sorted_histogram_of_biotypes(df_sorted_count_biotype, plot_title,
                                          1, plot_path + plot_file) #0: display
                                                                    #1: save file

    #
    # iterate throughout each gene type
    for gene_type in sortedByCounts_biotype:
        if 1 and (gene_type not in CHECK_GENE_TYPE): # just for several gene types
            continue

        if (gene_type not in list_of_all_biotypes):
            list_of_all_biotypes.append(gene_type)
            dict_of_all_df_species[gene_type] = df_species.copy()

        df_genes = df_all_genes[(df_all_genes["biotype"] == gene_type)]
        print("There are ", len(df_genes), " species of this type (", gene_type, ")", sep="")
        #
        # Calculate the difference of length to the closer gene
        #######################################################
        if 1:
            df_genes.reset_index(drop=True)
            pd.set_option('display.max_columns', None)
            l1 = df_genes['length'].to_list()
            l = np.array(l1)
            lprev = np.roll(l, 1)
            lnext = np.roll(l, -1)
            d = lnext-l # distance with the next
            d[-1] = 999999999
            dprev = l-lprev # distance with the prev
            dprev[0] = 999999999
            d[np.less(dprev,d)] = dprev[np.less(dprev,d)]
            df_genes["diffLength"] = d
            plot_path = OUTPUT_PROJECT_PATH+DATA_PATH_NAME+"/"
            plot_file_stem = DB +"_" + ENSEMBL_VERSION +"." + annotation_name + "." + gene_type
            df_genes["log10ofLength"] = np.log10(df_genes.length)
            if 0:
                pd.options.display.max_columns = None
                print(df_genes.head(5))
            benford_mod.save_columns_givenDf(df_genes, ["biotype", "length", "diffLength", "log10ofLength", "contig", "start", "end", "gene_id", "gene_name"],
                                             1, plot_path, plot_file_stem)
            del df_genes["diffLength"] # I do not know if this is necessary ...but just in case
            del d

        # Class:
        # DigitDistribution initialization
        ##################################
        gene_dist = first_digit_distributions.DigitDistribution(DB +":" + ENSEMBL_VERSION +":" + annotation_name + ":" + gene_type)
        gene_dist.fromLength__assig_count_and_prob(df_genes["length"].to_list(),
                                                       df_genes["gene_id"].to_list()) #1d and 2d
        if 0: #just for displaying the distribution
            gene_dist.print_DigitDistribution("all") #"1d", "2d" or "all"; "all" by default
        #
        # plot the length distribution for this gene type
        plot_title = DB.capitalize() +" " + ENSEMBL_VERSION +" (" + annotation_name + ")"
        x_axis_label = "length  ("+str(len(df_genes.length)) + " " + gene_type +")"
        plot_path = OUTPUT_PROJECT_PATH + GENE_TYPE_PATH_NAME + "/"+ gene_type + "/"
        plot_file = DB+"_"+ENSEMBL_VERSION+"."+annotation_name+"."+GENE_TYPE_PATH_NAME+"."+gene_type+".png"
        plot_dfColumn_histogram(df_genes, "length", plot_title, x_axis_label,
                              1, plot_path+plot_file)

        # initialize the Benford law distribution
        benford_dist = first_digit_distributions.DigitDistribution("benford law:::benford law")
        benford_dist.prob_of_benfordLaw("all")
        if 0:
            benford_dist.print_DigitDistribution("all") #1d, 2d, all

        # Comparisons
        #############
        stat_describe_length = df_genes[["length", "log10ofLength"]].describe()
        if 0:
            print(stat_describe_length)
            print(stat_describe_length["length"]["std"])
        dict_of_all_df_species[gene_type].loc[i, "count"]  = int(stat_describe_length["length"]["count"])
        #
        dict_of_all_df_species[gene_type].loc[i, "mean"]   = stat_describe_length["length"]["mean"]
        dict_of_all_df_species[gene_type].loc[i, "std"]    = stat_describe_length["length"]["std"]
        dict_of_all_df_species[gene_type].loc[i, "min"]    = int(stat_describe_length["length"]["min"])
        dict_of_all_df_species[gene_type].loc[i, "25perc"] = stat_describe_length["length"]["25%"]
        dict_of_all_df_species[gene_type].loc[i, "50perc"] = stat_describe_length["length"]["50%"]
        dict_of_all_df_species[gene_type].loc[i, "75perc"] = stat_describe_length["length"]["75%"]
        dict_of_all_df_species[gene_type].loc[i, "max"]    = int(stat_describe_length["length"]["max"])
        #
        dict_of_all_df_species[gene_type].loc[i, "log10_mean"]   = stat_describe_length["log10ofLength"]["mean"]
        dict_of_all_df_species[gene_type].loc[i, "log10_std"]    = stat_describe_length["log10ofLength"]["std"]
        dict_of_all_df_species[gene_type].loc[i, "log10_min"]    = stat_describe_length["log10ofLength"]["min"]
        dict_of_all_df_species[gene_type].loc[i, "log10_25perc"] = stat_describe_length["log10ofLength"]["25%"]
        dict_of_all_df_species[gene_type].loc[i, "log10_50perc"] = stat_describe_length["log10ofLength"]["50%"]
        dict_of_all_df_species[gene_type].loc[i, "log10_75perc"] = stat_describe_length["log10ofLength"]["75%"]
        dict_of_all_df_species[gene_type].loc[i, "log10_max"]    = stat_describe_length["log10ofLength"]["max"]
        #
        dfD_1d, dfB_1d = prepare_plot_compare2Benford(gene_dist, benford_dist, "1d")
        dfD_2d, dfB_2d = prepare_plot_compare2Benford(gene_dist, benford_dist, "2d")
        plot_path = OUTPUT_PROJECT_PATH+COMPARE_PATH_NAME+"/"
        #
        # chiSqrt
        #########
        if 1:
            plot_file_stem = "/"+DB+"_"+ENSEMBL_VERSION+"."+annotation_name+"."+gene_type+"."+CHI_SQRT_PATH_NAME
            #
            x_sqrt, pvalue, dof = plot_singleComparison2Benford(dfD_1d, dfB_1d, gene_dist.name,
                                                                "1d", "point-line", "count", 1, 1,  # "1d" or "2d"
                                                                plot_path + CHI_SQRT_PATH_NAME + plot_file_stem +
                                                                "_1d" + ".png")  # "point-line" or "bar"
                                                                                # "count" or "prob" )
                                                                                # bool_plot_line: 0 or 1
                                                                                # bool_plot: 0 or 1
            dict_of_all_df_species[gene_type].loc[i, "1d_xs_xsqrt"] = x_sqrt
            dict_of_all_df_species[gene_type].loc[i, "1d_xs_pvalue"] = pvalue
            dict_of_all_df_species[gene_type].loc[i, "1d_xs_dof"] = dof
            #
            x_sqrt, pvalue, dof = plot_singleComparison2Benford(dfD_2d, dfB_2d, gene_dist.name,
                                                "2d", "point-line", "count", 1, 1,
                                                plot_path + CHI_SQRT_PATH_NAME + plot_file_stem +
                                                "_2d" +".png") # "1d" or "2d"
                                                              # "point-line" or "bar"
                                                              # "count" or "prob" )
                                                              #  bool_plot_line: 0 or 1
                                                              #  bool_plot: 0 or 1
            dict_of_all_df_species[gene_type].loc[i, "2d_xs_xsqrt"] = x_sqrt
            dict_of_all_df_species[gene_type].loc[i, "2d_xs_pvalue"] = pvalue
            dict_of_all_df_species[gene_type].loc[i, "2d_xs_dof"] = dof
        #
        # correlation
        #############
        if 1:
            plot_file_stem = "/" + DB +"_" + ENSEMBL_VERSION +"." + annotation_name +"." + gene_type +"." + CORR_PATH_NAME
            pearson, pval = plot_2var_correlation(dfD_1d, dfB_1d, gene_dist.name, "1d",                   # "1d" or "2d"
                                                  1, plot_path+CORR_PATH_NAME+plot_file_stem+"_1d"+".png")# bool_plot_line: 0 or 1
                                                                                                          # bool_plot: 0 or 1
            dict_of_all_df_species[gene_type].loc[i, "1d_corr_pearson"] = pearson
            dict_of_all_df_species[gene_type].loc[i, "1d_corr_pvalue"] = pval
            pearson, pval = plot_2var_correlation(dfD_2d, dfB_2d, gene_dist.name, "2d",                   # "1d" or "2d"
                                                  1, plot_path+CORR_PATH_NAME+plot_file_stem+"_2d"+".png")# bool_plot_line: 0 or 1
                                                                                                          # bool_plot: 0 or 1
            dict_of_all_df_species[gene_type].loc[i, "2d_corr_pearson"] = pearson
            dict_of_all_df_species[gene_type].loc[i, "2d_corr_pvalue"] = pval
        #
        # Nigrini
        #############
        m_ad, max_ad = print_results_Nigrini2000(dfD_1d["prob"].to_list(),
                                  dfB_1d["prob"].to_list(), 9, 1) # 9 dims and 1:print
        dict_of_all_df_species[gene_type].loc[i, "1d_nigrini2000Mad", ] = m_ad
        dict_of_all_df_species[gene_type].loc[i, "1d_nigrini2000_MAXad"] = max_ad
        if 0:
            m_ad, max_ad = print_results_Nigrini2000(dfD_2d["prob"].to_list(),
                                        dfB_2d["prob"].to_list(), 90, 1) # 9 dims and 1:print

plot_file_grund = OUTPUT_PROJECT_PATH+DATA_PATH_NAME+"/"
for b in list_of_all_biotypes:
    plot_file = plot_file_grund+DB+"_"+ENSEMBL_VERSION+".species."+DATA_PATH_NAME+"."+b+".tsv"
    dict_of_all_df_species[b].to_csv(plot_file, index=False, sep="\t") #log, only the first indexing

if BOOL_PLOT_GENETYPE_HIST:
    geneTypePerSpecie_file = plot_file_grund + DB + "_" + ENSEMBL_VERSION + "." + GENE_TYPE_PER_SPECIE + ".tsv"
    df_biotypes_per_specie.to_csv(geneTypePerSpecie_file, index=False, sep="\t", na_rep='None')
if 0:
    print(df_biotypes_per_specie)