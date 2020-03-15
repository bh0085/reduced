import os, re, sys
import pandas as pd
import datetime

PRJ_DIR = '/data/cgs/bh0085/prj/reduced_all/reduced_1219'
PRJ_NAME = 'reduced_1219'

today = datetime.date.today()
DATESTR = "{0:02}{0:02}".format(today.month,today.day)

#######################################################
# Note: Directories should end in / always
#######################################################
SRC_DIR = os.path.join(PRJ_DIR , 'src/')
DATA_DIR = os.path.join(PRJ_DIR , 'data/')
LIBS_DIR = os.path.join(PRJ_DIR , 'libs/')
OUT_PLACE = os.path.join(PRJ_DIR , 'out/')
QSUBS_DIR = os.path.join( PRJ_DIR , 'qsubs/')
FIGS_PLACE = os.path.join(PRJ_DIR , 'figs/')
LOGS_DIR = os.path.join(PRJ_DIR, 'logs/')


#dictionary keyed by sequencing group names, containing:
#  "reads_place" (with folders for each sequencing lane)
#  "lib_file" (having names and descriptions associated with each sequencing lane)
SEQUENCING_READS_META = {
                 "reads_1219":{"reads_place":'/data/cgs/bh0085/shortreads/1219_SEQUENCING_READS',
                            "lib_file":os.path.join(LIBS_DIR,"1bp_SaCasgRNA_redlibrary_SequenceSamples.csv")}
               }


reads_dir = "/data/cgs/bh0085/shortreads/1219_SEQUENCING_READS/"
dm_reads_dir = "/data/cgs/bh0085/shortreads/1219_SEQUENCING_READS/READS_DEMULTIPLEXED/"

if not os.path.isdir(dm_reads_dir):
    os.makedirs(dm_reads_dir)
    
exp_lib = pd.read_csv(os.path.join("../libs/1bp_SaCasgRNA_redlibrary_SequenceSamples.csv"))
exp_idx_code = exp_lib.apply(lambda x:f"{x['Index 1 (i7) sequence']}:{x['Index 2 (i5) sequence']}",axis=1)
exp_lib["idx_code"] = exp_idx_code   
exp_lib = exp_lib.set_index("idx_code").sort_index()
exp_lib["exptype"] =exp_lib.apply(lambda x: "+ku60019" if "+" in x["Sample name"] else "-ku60019" ,axis=1)

exp_lib["r1_file"] = exp_lib.apply(lambda x:os.path.join(dm_reads_dir,f"{x.exptype}_{x.rep}_R1.fastq"),axis=1)
exp_lib["r2_file"] = exp_lib.apply(lambda x:os.path.join(dm_reads_dir,f"{x.exptype}_{x.rep}_R2.fastq"),axis=1)


reads_place = dm_reads_dir
reads_files = os.listdir(reads_place)
info = exp_lib

info["r1_fastq_path"] = info.r1_file
info["r2_fastq_path"] = info.r2_file

info["nm"] = "all"
info["reads_place"] = reads_place
info["Name"] =  info["Sample name"]
info["Description"] = info["exptype"]

SEQUENCING_INFO = info     
SEQUENCING_INFO["replicate"] = SEQUENCING_INFO["rep"]


LIBRARY_DF = pd.read_csv("/data/cgs/bh0085/prj/reduced_all/reduced_1219/libs/1219_sacas9_library.csv")
print(LIBRARY_DF["Identifier number"])
LIBRARY_DF["Name"] = LIBRARY_DF["Identifier number"].apply(lambda x: int(x)) - 1

LIBRARY_DF.set_index("Name")

REDUCED_LIB = LIBRARY_DF
