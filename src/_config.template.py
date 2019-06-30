import os, re, sys
import pandas as pd
import datetime

PRJ_DIR = '/cluster/bh0085/prj/reduced_all/heysol'
PRJ_NAME = 'REDUCEDLIB_heysol_June2019'

today = datetime.date.today()
DATESTR = "{0:02}{0:02}".format(today.month,today.day)

#######################################################
# Note: Directories should end in / always
#######################################################
SRC_DIR = os.path.join(PRJ_DIR , 'src/')
DATA_DIR = os.path.join(PRJ_DIR , 'data/')
OUT_PLACE = os.path.join(PRJ_DIR , 'out/')
QSUBS_DIR = os.path.join( PRJ_DIR + 'qsubs/')
FIGS_PLACE = os.path.join(PRJ_DIR , 'figs/')
LOGS_DIR = os.path.join(PRJ_DIR, 'logs/')


#parse sequencing data information 
from a0_library_loaders_heysol_subproj import SEQUENCING_INFO, LIBRARY_DF
