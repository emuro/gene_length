### import sys
### sys.path.append('./lib/')
from lib import constants as c
from lib import os_utils

from pyensembl import Genome
import pandas as pd
import re  # for regular expressions
import glob, os  # to deal with the dirs and files
import time

""" 
    IMPORTANT NOTE:
        -install pyensembl
            It was easy to install it in Pycharm. But for downloading the data I did the next:
            (see https://pypi.org/project/pyensembl/)
            * pip install pyensembl
            * pyensembl install --release 98 --species human
"""


def indexData_for_all_species(ensembl_division="vertebrates", ensembl_version=49,
                              db_local_path_root="/Users/enriquem.muro/data/compressed/",
                              path="/home/emuro/data/compressed/ftp.ensembl.org/pub/metazoa/release-49/gtf"):
    ''' given the next inputs:
            -ensembl_version: (ie. 98)
            - db_path: (where is the ensembl local data?)

        It index all the gene annotation for all the species:
            gtf gene annotations:
                there is a directory for each annotated species and each one
                (ie. Homo_sapiens) contains a compressed file (*gtf.gz)
                with the ensembl gene annotations

                (ie. Homo_sapiens.GRCh38.98.gtf.gz). The file name has the
                next format: "species"."assembly"."ensembl_version"."gtf.gz"

        and then it performs the indexing It created a file *db within each species dir

         output:
     '''
    species=[]
    species_withpath=[]
    os_utils.recursiveSearch_leave_subdirs(path, species_withpath)  # output in path_list
    for p in species_withpath:
        species.append(os.path.basename(p))
    if 0:
        print(species)
        print("there are", len(species), "species annotated for", ensembl_division)

    endOfFile_in_unix = "." + str(ensembl_version) + ".gtf.gz"
    df_species = pd.DataFrame([],
                              columns=["species", "assembly", "ensembl_division", "ensembl_version", "path", "db_file", "file"])

    for s, s_path in zip(species, species_withpath):
        os.chdir(s_path)
        for file in glob.glob("*" + endOfFile_in_unix):
            if file in c.GTF_FILES_TO_EXCLUDE:
                continue
            annot_file = s.capitalize() + ".(.*)" + endOfFile_in_unix #reg exp
            res = re.findall(annot_file, file)
            # finally index the data base; check the time it takes
            start = time.time()
            new_species_row = {"species":s, "assembly":res[0],
                               "ensembl_division":ensembl_division, "ensembl_version":ensembl_version,
                               "db_path":re.sub(r"{}".format(db_local_path_root), '', s_path),
                               "db_file":re.sub(r'.gtf.gz', '.gtf.db', file),
                               "file":file,
                               "indexing_time":-999.666}
            if 0:
                print(new_species_row)
                print(file)
            index_db_of_geneAnnotation_from_gtfFile(new_species_row["assembly"],
                                                    new_species_row["species"],
                                                    new_species_row["file"])
            new_species_row["indexing_time"]=time.time()-start
            df_species=df_species.append(new_species_row, ignore_index=True)
    return df_species



def indexData_for_all_species__all_gz_in_a_dir(gz_path="/Users/enriquem.muro/data/compressed/geneLength/ensemblgenomes_98_and_49/"):
    ''' given the next inputs:
            - gz_path: (where is the ensembl local data? (all the files *gz that I previously)

        It index all the gene annotation for all the species:
            gtf gene annotations:
                there is a file for each annotated species and each one
                (ie. Homo_sapiens) contains a compressed file (*gtf.gz)
                with the ensembl gene annotations

                (ie. Homo_sapiens.GRCh38.98.gtf.gz). The file name has the
                next format: "species"."assembly"."ensembl_version"."gtf.gz"

        and then it performs the indexing It created a file *db within each species dir

        output:
        a data frame with the information of the files that have been indexed
     '''
    #
    # get the species
    os.chdir(gz_path)
    fileExtension=".gtf.gz"
    print(gz_path)
    species=glob.glob("*"+fileExtension, recursive=False)
    species.sort(key=os.path.getmtime) #sort the files by modidification time
    #endOfFile is ".gtf.gz": take care about which is it, no chr, no abinitio, etc.
    df_species=pd.DataFrame([],columns=["species","assembly","ensembl_version","path","file"])
    for file in species:
        if 0:
            print (file)

        #dirty trick to get species, genome, version...:-(
        aux=file.split(".gtf.gz")[0]
        aux2=aux.split(".")
        v=aux2.pop() #ensembl version
        s=aux2[0] #species
        del aux2[0]
        a=".".join(aux2) #assembly

        # finally index the data base; check the time it takes
        start=time.time()
        new_species_row={"species": s,"assembly": a,"ensembl_version": v,
                         "path": gz_path,"file": file,
                         "indexing_time": -999.666}
        data=index_db_of_geneAnnotation_from_gtfFile(new_species_row["assembly"],
                                                     new_species_row["species"],
                                                     new_species_row["file"])
        new_species_row["indexing_time"]=time.time()-start
        df_species=df_species.append(new_species_row,ignore_index=True)
    return df_species



def index_db_of_geneAnnotation_from_gtfFile(ref_name="GRCh38", annot_name="homo_sapiens",
        gtf_file="/home/emuro/data/compressed/ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.gtf.gz"):
    ''' given the next inputs:
            ref_name: the reference version, the genome version
            annot_name: the species
            gtf file: the annotation file (compressed *gz or not) with path
        index the db annotation file
        output:
            data: with the indexed db object
    '''

    data = Genome(  # parse GTF and construct database of genomic features
        reference_name  = ref_name,
        annotation_name = annot_name,
        gtf_path_or_url = gtf_file)
    data.index() # parse GTF and construct database of genomic features
    return data



def retrieve_geneAnnotationObject_fromIndexedData(data):
    ''' given the next inputs:
           data: already indexed data
        retrieve the gene annotation of all the biotypes
        output:
           gene object
    '''

    genes = data.genes(contig=None, strand=None)
    return genes



def geneObjectList_to_df_sorted_by_biotype_and_geneLength(genes):
    ''' given the next inputs:
            -species:
            list of gene class instances obtained previously with
            get_all_genes_of_biotype. That is, the list of gene
            annotations for a gene type (aka biotype)
        output:
            df_all_genes (a data frame with the gene annotation of all the species)
            sorted by: first by biotype and then by gene length.
            It deletes columns: db and genome
    '''

    g=genes[0] # get the diff annotations of the gene
    attributes = [a for a, v in g.__dict__.items()
            if not re.match('<function.*?>', str(v))
            and not (a.startswith('__') and a.endswith('__'))]

    dict_genes=dict()
    for att in attributes: #dict initialization
        dict_genes[att] = []
    for g in genes: #annotation of each species of the biotype
        for att in attributes:
            dict_genes[att].append(getattr(g, att))  #update dict
    df_all_genes=pd.DataFrame(dict_genes)
    # deleting columns: genome and df. Not necessary
    del df_all_genes["genome"]
    del df_all_genes["db"]
    df_all_genes["length"] = df_all_genes["end"] - df_all_genes["start"]
    df_all_genes = df_all_genes.sort_values(by=["biotype", "length"])
    if 0: # just show the df (all species captured by biotype)
        pd.set_option('display.max_columns', None)
        print(df_all_genes)
        print(df_all_genes.dtypes)
    if 0:
        print("\tAfter use_pyensembl.geneObjectList_to_df_sorted_by_biotype_and_geneLength, ", "there are", len(df_all_genes),
              "genes (df_all_genes)")
    return df_all_genes



def count_genes_by_biotype_givenDf(df_all_genes,
                                   min_number_of_genes_per_biotype=50):
    ''' given the next inputs:
            ENSEMBL_VERSION
            #SPECIES [Human at the moment]
            MIN_NUMBER_OF_GENES_PER_BIOTYPE (int)
        retrieve a df with the number of species for each gene type (aka biotype)
        output:
            df_count_of_each_biotype #sorted by the counts of each biotype
    '''
    if 0:
        print(len(df_all_genes), "number of species")
        print(df_all_genes.head(10))
        exit()
    # check how many species are annotated for each type of gene (aka biotype)
    #
    numOfGenes_per_biotype={}
    numOfGenes_per_contig={}
    for index, row_gene in df_all_genes.iterrows():
        biotype = row_gene["biotype"]
        if biotype in numOfGenes_per_biotype.keys():
            numOfGenes_per_biotype[biotype] += 1
        else:
            numOfGenes_per_biotype[biotype] = 1
        if 0:
            contig = row_gene["contig"]
            if contig in numOfGenes_per_contig.keys():
                numOfGenes_per_contig[contig] += 1
            else:
                numOfGenes_per_contig[contig] = 1
    if 0:
        print(len(numOfGenes_per_biotype), "diff. types of species (aka biotypes)")
        print(numOfGenes_per_biotype)
    if 0:
        print(len(numOfGenes_per_contig), "diff. types of chrs (aka contig)")
        print(numOfGenes_per_contig)
        print(sorted(numOfGenes_per_contig.keys()))

    sorted_numOfGenes_per_biotype = sorted(numOfGenes_per_biotype.items(),   #sort by counts
                                           key=lambda x: x[1], reverse=True) #list, descending
    if 0:
        print(sorted_numOfGenes_per_biotype) # list of list
    if 0:
        sorted_numOfGenes_per_contig = sorted(numOfGenes_per_contig.items(),
                                               key=lambda x: x[1], reverse=True) #list
        print(sorted_numOfGenes_per_contig) # list of list
    #
    # prepare a df with the counts (of biotype)
    data=dict
    aux1=list()
    aux2=list()
    for b in sorted_numOfGenes_per_biotype:
        if (b[1]>min_number_of_genes_per_biotype):
            aux1.append(b[0])
            aux2.append(b[1])
    data = {'biotype':  aux1,
            'number': aux2
            }
    df = pd.DataFrame (data, columns = ['biotype', 'number'])
    if 0:
        print(df)
    return df
