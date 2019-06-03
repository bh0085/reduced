import os, re, sys
import pandas as pd
import datetime

PRJ_DIR = '/cluster/bh0085/prj/reduced_all/reduced_sm/'
PRJ_NAME = 'REDUCEDLIB_SMALLMOL_May2019'

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


#dictionary keyed by sequencing group names, containing:
#  "reads_place" (with folders for each sequencing lane)
#  "lib_file" (having names and descriptions associated with each sequencing lane)
SEQUENCING_READS_META = {
    #"reads_28": {"reads_place":"/cluster/bh0085/shortreads/190328Gif",
                #              "lib_file":"../libs/smallmol2 sequencing directory contents - lib28.csv"},
                # "reads_29":  {"reads_place":"/cluster/bh0085/shortreads/190329Gif",
                #               "lib_file":"../libs/smallmol2 sequencing directory contents - lib29.csv"},
                 "reads_502":{"reads_place":'/cluster/bh0085/shortreads/190502Gif',
                            "lib_file":"../libs/smallmol2 sequencing directory contents - lib502.csv"}
               }


SEQUENCING_INFO = pd.DataFrame()
for nm,vals in SEQUENCING_READS_META.items():
    #read sequencing maps and experimental design files
    reads_place = vals["reads_place"]
    lib_file = vals["lib_file"]
    reads_files = os.listdir(reads_place)
    info = pd.read_csv(lib_file) 

    #add metadata to sequencing info DF
    info = pd.concat([info,
                      info.apply(
        lambda x:pd.Series( re.compile("(?P<Name>D[^\s]+) : (?P<Description>[^()]*)\((?P<Index>[^)]*)\)")
                           .search(x["code"]).groupdict()),axis=1)]
                     ,axis=1)
    info["nm"] = nm
    info["reads_place"] = reads_place
    reads_files =os.listdir(reads_place)
    
    info["reads_location"] = info.Name.apply(
        lambda x: [e for e in reads_files if x in e]).apply(
            lambda x:os.path.join(reads_place,x[0]) if len(x)>0 
        else None)
    
    #find fastq files using sequencing info DF
    if False in set(pd.notna(info.reads_location)):
        raise Exception("missing reads data file for {0}".format(nm))
    info["fastq"] = info.apply(lambda x:[fn for fn in os.listdir(x.reads_location) if ~(pd.isna(x.reads_location)) & ( fn[-5:]=="fastq") & ("_2_" in fn)][0],axis=1)
    info["fastq_path"] = info["reads_location"] +"/"+ info["fastq"]
    SEQUENCING_INFO = SEQUENCING_INFO.append(info)
     
 
SEQUENCING_INFO["Description"] = SEQUENCING_INFO.apply(lambda x:re.compile(": (.*)\(.*\)").search(x.code).groups()[0].strip(),axis=1)
SEQUENCING_INFO = SEQUENCING_INFO.loc[18:77]
SEQUENCING_INFO["drug_name"] = SEQUENCING_INFO.Description.str.slice(0,-2)
SEQUENCING_INFO["replicate"] = SEQUENCING_INFO.Description.str.slice(-1)

REDUCED_LIB = pd.read_csv("/cluster/bh0085/data/reduced_lib_design.csv")
REDUCED_LIB["Name"] = REDUCED_LIB["Identifier number"].apply(lambda x: int(x)) - 1
REDUCED_LIB.set_index("Name")

