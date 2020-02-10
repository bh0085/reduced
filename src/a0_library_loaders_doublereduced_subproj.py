import pandas as pd
import os


#dictionary keyed by sequencing group names, containing:
#  "reads_place" (with folders for each sequencing lane)
#  "lib_file" (having names and descriptions associated with each sequencing lane)
SEQUENCING_READS_META = {
                "doublereduced":{"reads_place":"/cluster/bh0085/shortreads/FC_04530/Unaligned_1234_PF_mm1/Data/Project_richsherwood/",
                "lib_file":"../libs/FC04530_doublereduced_readsmap.csv"
                                }
               }

SEQUENCING_INFO = pd.DataFrame()
for nm,vals in SEQUENCING_READS_META.items():
    reads_place = vals["reads_place"]
    lib_file = vals["lib_file"]


    reads_files = os.listdir(reads_place)
    info = pd.read_csv(lib_file) 

#     #add metadata to sequencing info DF
#     info = pd.concat([info,
#                     info.apply(
#         lambda x:pd.Series( re.compile("(?P<Name>D[^\s]+) : (?P<Description>[^()]*)\((?P<Index>[^)]*)\)")
#                         .search(x["code"]).groupdict()),axis=1)]
#                     ,axis=1)
    info["code_offset"] = info.code.apply(lambda x: int(x[-3:]))
    info["s_num"] = info["code_offset"] - 688
    lanes = pd.DataFrame([pd.Series({"lane_num":i}) for i in range(1,5)])
    info = pd.merge(info.assign(key=1),lanes.assign(key=1),on="key")
    info["r1_fastq_path"] = info.apply(lambda x: os.path.join(reads_place,f"LIB037937_{x.code}_S{x.s_num}_{x.lane_num}_R1.fastq"),axis=1)
    info["r2_fastq_path"] = info.apply(lambda x: os.path.join(reads_place,f"LIB037937_{x.code}_S{x.s_num}_{x.lane_num}_R1.fastq"),axis=1)
    info["nm"] = nm
    info["reads_place"] = reads_place
    info["Name"] =  info.apply(lambda x:f"LIB037937_{x['code']}_S{x['s_num']}_L00{x['lane_num']}_R2",axis = 1)
    info["Description"] = info["desc"]

    SEQUENCING_INFO = SEQUENCING_INFO.append(info)


#LOAD THE LIBRARY CSV
REDUCED_LIB = pd.read_csv("/cluster/bh0085/data/reduced_lib_design.csv")
REDUCED_LIB["Name"] = REDUCED_LIB["Identifier number"].apply(lambda x: int(x)) - 1
REDUCED_LIB.set_index("Name")

BIGGER_LIB = pd.merge(REDUCED_LIB[["Name"]].rename({"Name":"Name1"},axis="columns").assign(key=1),REDUCED_LIB[["Name"]].rename({"Name":"Name2"},axis="columns").assign(key=1))
BIGGER_LIB["seqleft"] =BIGGER_LIB.Name1.apply(lambda x: REDUCED_LIB.loc[x].seqleft)
BIGGER_LIB["seqright"] =BIGGER_LIB.Name2.apply(lambda x: REDUCED_LIB.loc[x].seqright)
BIGGER_LIB["targetseq_61"] = BIGGER_LIB["seqleft"] + BIGGER_LIB["seqright"]
BIGGER_LIB["StringName"] = BIGGER_LIB.apply(lambda x: f"{x.Name1}-{x.Name2}",axis = 1)
BIGGER_LIB["Name"] = BIGGER_LIB.apply(lambda x: f"{int(x.Name1)*48+int(x.Name2)}",axis = 1)
BIGGER_LIB = BIGGER_LIB.drop("key",axis = "columns")

LIBRARY_DF= BIGGER_LIB