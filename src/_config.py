import sys
import os

PRJ_DIR = '/cluster/bh0085/prj/reduced2/'
SRC_DIR = PRJ_DIR + 'src/'


#######################################################
# Note: Directories should end in / always
#######################################################
DATA_DIR = PRJ_DIR + 'data/'
OUT_PLACE = PRJ_DIR + 'out/'
RESULTS_PLACE = PRJ_DIR + 'results/'
QSUBS_DIR = PRJ_DIR + 'qsubs/'
READS_DIR = '/cluster/bh0085/shortreads/FC_04530/Unaligned_1234_PF_mm1/Data/Project_richsherwood'

EXP_DESIGN_FILE = os.path.join(DATA_DIR, "reduced.csv")

#######################################################
#######################################################

CLEAN = False       # Values = 'ask', True, False

