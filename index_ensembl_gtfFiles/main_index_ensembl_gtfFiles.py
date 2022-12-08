#
# DB INDEXING
# NOTE: NEEDS TO BE RUN ONLY ONE TIME. That is, only one indexing per *.gtf
# pyemsembl:
# Index the gene annotation file (*.gtf) for the species given an ensembl version.
# A *.gtf file is indexed in ~30sec per species. A file needs to be indexed only once time.

import sys
sys.path.append('./lib/')
from lib import constants as c  # (or constants__arcturus.py)
from lib import use_pyensembl as pyens

if 1: # index the encode annotation files for all the species
    BOOL_FIRST_INDEXING = 1
    df_species = pyens.indexData_for_all_species(c.ENSEMBL_DIVISION, c.ENSEMBL_VERSION,
                                                 c.DB_LOCAL_PATH_ROOT, c.GTF_LOCAL_PATH)
    if BOOL_FIRST_INDEXING:
        df_species.to_csv(c.LOG_ANNOTATION_INDEXED_SPECIES_FILE, index = False, sep="\t") #log, only the first indexing



# #
# # DB INDEXING
# # NOTE: NEEDS TO BE RUN ONCE PER SESSION
# # pyemsembl: Index the gene annotation for all the species given an
# #            ensembl version. Once that it does it (~30 sec per species)
# #            in the next session is extremely fast (secs for all species)
# if 0: # index the encode annotation files for all the species
#     BOOL_FIRST_INDEXING = 0
#     df_species=use_pyensembl.indexData_for_all_species(ENSEMBL_VERSION, DB_LOCAL_PATH)
#     if BOOL_FIRST_INDEXING:
#         df_species.to_csv(LOG_ANNOTATION_INDEXED_SPECIES_FILE, index = False, sep="\t") #log, only the first indexing
#     exit()


