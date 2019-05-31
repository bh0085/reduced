import os, re, sys
import pandas as pd


PRJ_DIR = '/cluster/bh0085/prj/reduced/reduced_sm/'


#######################################################
# Note: Directories should end in / always
#######################################################
SRC_DIR = os.path.join(PRJ_DIR , 'src/')
DATA_DIR = os.path.join(PRJ_DIR , 'data/')
OUT_PLACE = os.path.join(PRJ_DIR , 'out/')
QSUBS_DIR = os.path.join( PRJ_DIR + 'qsubs/')


#dictionary keyed by sequencing group names, containing:
#  "reads_place" (with folders for each sequencing lane)
#  "lib_file" (having names and descriptions associated with each sequencing lane)
SEQUENCING_READS_META = {"reads_28": {"reads_place":"/cluster/bh0085/shortreads/190328Gif",
                             "lib_file":"../libs/smallmol2 sequencing directory contents - lib28.csv"},
                "reads_29":  {"reads_place":"/cluster/bh0085/shortreads/190329Gif",
                              "lib_file":"../libs/smallmol2 sequencing directory contents - lib29.csv"},
                "reads_502":{"reads_place":'/cluster/bh0085/shortreads/190502Gif',
                            "lib_file":"../libs/smallmol2 sequencing directory contents - lib502.csv"}
               }

REDUCED_LIB = pd.read_csv("~/data/reduced_lib_design.csv")
REDUCED_LIB["Name"] = REDUCED_LIB["Identifier number"].apply(lambda x: int(x)) - 1
REDUCED_LIB.set_index("Name")

