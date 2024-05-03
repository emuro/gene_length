# python3
# ################################################################## #
# retrieve_data__functions.py (C) 2023 Mainz.
# Author: Enrique M. Muro
# ################################################################## #
#
# --------------------------------------------------------------------
# Project: gene_length
#
# ################################################################## #


import numpy as np
import pandas as pd

import sys
sys.path.append('./')
import gl_constants as gl_c



def retrieve_species_with_high_quality_genomes(): # HQG: high quality genome
    highQ_df = pd.read_csv(gl_c.G_NCBI_GENOME_DATA_FILE,
                           low_memory=False, sep="\t") # tax_id, status, accession  
    highQ_df = highQ_df[["species","assembly_status"]]
    if 0:
        display(highQ_df.head(2))
        display(highQ_df.columns)
    
    # status can be:["Complete Genome", "Chromosome", "Scaffold", "Contig"]
    good_status = ["Complete Genome", "Chromosome"] 

    # filter only to high Quality genomes
    highQ_df = highQ_df[highQ_df['assembly_status'].isin(good_status)]
    display(highQ_df.shape)
    return highQ_df["species"].tolist()


def from_species_list_retrieve_species_with_HQG(species_l): # HQG: high quality genome
    species_s = set(species_l)
    hqg_s     = set(retrieve_species_with_high_quality_genomes())
    return list(species_s.intersection(hqg_s))
