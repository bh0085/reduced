import pandas as pd
import os



#dictionary keyed by sequencing group names, containing:
#  "reads_place" (with folders for each sequencing lane)
#  "lib_file" (having names and descriptions associated with each sequencing lane)
SEQUENCING_READS_META = {
                "heysol":{"reads_place":"/cluster/bh0085/shortreads/FC_05081/Unaligned_1234_PF_mm1/Data/Project_richsherwood",
                "lib_file":os.path.join("../libs/0630_reducedlib_heysol_sequencing_maps.csv")
                                }
               }

SEQUENCING_INFO = pd.DataFrame()
for nm,vals in SEQUENCING_READS_META.items():
    reads_place = vals["reads_place"]
    lib_file = vals["lib_file"]


    reads_files = os.listdir(reads_place)
    info = pd.read_csv(lib_file) 

    info["code_offset"] = info.code.apply(lambda x: int(x[-3:]))
    info["s_num"] = info["code_offset"] - 83
    lanes = pd.DataFrame([pd.Series({"lane_num":i}) for i in range(1,5)])
    info = pd.merge(info.assign(key=1),lanes.assign(key=1),on="key")
    info["r1_fastq_path"] = info.apply(lambda x: os.path.join(reads_place,f"LIB041331_{x.code}_S{x.s_num}_L00{x.lane_num}_R1.fastq"),axis=1)
    info["r2_fastq_path"] = info.apply(lambda x: os.path.join(reads_place,f"LIB041331_{x.code}_S{x.s_num}_L00{x.lane_num}_R2.fastq"),axis=1)
    info["nm"] = nm
    info["reads_place"] = reads_place
    info["Name"] =  info.apply(lambda x:f"LIB041331_{x['code']}_S{x['s_num']}_L00{x['lane_num']}",axis = 1)
    info["Description"] = info["desc"]

    SEQUENCING_INFO = SEQUENCING_INFO.append(info)


#LOAD THE LIBRARY CSV
REDUCED_LIB = pd.read_csv("/cluster/bh0085/data/reduced_lib_design.csv")
REDUCED_LIB["Name"] = REDUCED_LIB["Identifier number"].apply(lambda x: int(x)) - 1
REDUCED_LIB.set_index("Name")
LIBRARY_DF= REDUCED_LIB