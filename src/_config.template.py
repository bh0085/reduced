import os, re, sys
import pandas as pd
import datetime

PRJ_DIR = '/cluster/bh0085/prj/reduced_all/doublereduced2'
PRJ_NAME = 'REDUCEDLIB_DOUBLEREDUCED_June2019'

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
READS_DIR="/cluster/bh0085/shortreads/FC_04530/Unaligned_1234_PF_mm1/Data/Project_richsherwood/"


#parse sequencing data information 
from a0_library_loaders_doublereduced_subproj import SEQUENCING_INFO, LIBRARY_DF
