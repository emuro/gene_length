# python3
# ################################################################## #
# main_rest_ensembl.py (C) Jan-2021 Mainz.
# Author: Enrique M. Muro
# ################################################################## #
#
# ------------------------------------------------------------------------
# Project:
#
# Use the next modules:
#   os_utils.py
#
#
import os

import sys
sys.path.append('./lib/')
from lib import os_utils

# search leave' subdirs from root
#################################
BOOL_SEARCH_LEAVE_SUBDIRS = 1
if BOOL_SEARCH_LEAVE_SUBDIRS:
    current_path="/Users/enriquem.muro/data/compressed/ftp.ensemblgenomes.org/pub/metazoa/release-49/gtf"
    path_leave_dirs=[]
    os_utils.recursiveSearch_leave_subdirs(current_path, path_leave_dirs) #output in path_list
    print("_________")
    #
    leave_dirs=[]
    for p in path_leave_dirs:
        leave_dirs.append(os.path.basename(p))
    print(leave_dirs)
    print("there are",len(leave_dirs),"leave dirs")

BOOL_UNIX_BASENAME = 0
if BOOL_UNIX_BASENAME:
    path_dir="/Users/enriquem.muro/data"
    paths=["/Users/enriquem.muro","/Users/enriquem.muro/data"]
    out_dirs=[]
    for p in paths:
        out_dirs.append(os.path.basename(p))
    print(out_dirs)