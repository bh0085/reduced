# 
from _config import SEQUENCING_READS_META, REDUCED_LIB, OUT_PLACE


import pandas as pd
import numpy as np
import os, sys, re

import sys, os, datetime, subprocess, imp
sys.path.append('/cluster/bh0085/')
from mybio import util


# Default params
NAME = util.get_fn(__file__)
out_dir = _config.OUT_PLACE + NAME + '/'
util.ensure_dir_exists(out_dir)


sequencing_info = pd.DataFrame()
for nm,vals in sequencing_reads_meta.items():
    #read sequencing maps and experimental design files
    reads_place = vals["reads_place"]
    lib_file = vals["lib_file"]
    reads_files = os.listdir(reads_place)
    #lib = pd.read_csv(lib_file,names=["Code"])
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
        raise Exception(f"missing reads data file for {nm}")
    info["fastq"] = info.apply(lambda x:[fn for fn in os.listdir(x.reads_location) if ~(pd.isna(x.reads_location)) & ( fn[-5:]=="fastq")][0],axis=1)
    info["fastq_path"] = info["reads_location"] +"/"+ info["fastq"]
    sequencing_info = sequencing_info.append(info)
     
 

##
# Functions
##
def split(inp_fn, out_nm):
  #print inp_fn
  inp_fn_numlines = util.line_count(inp_fn)

  #print inp_fn_numlines
  num_splits = 15
  split_size = int(inp_fn_numlines / num_splits)
  if num_splits * split_size < inp_fn_numlines:
    split_size += 1
  while split_size % 4 != 0:
    split_size += 1
  print (f'Using split size {split_size}')

  split_num = 0
  for idx in range(1, inp_fn_numlines, split_size):
    start = idx
    end = start + split_size  
    out_fn = out_dir + out_nm + '_%s.fastq' % (split_num)
    command = 'tail -n +%s %s | head -n %s > %s' % (start, inp_fn, end - start, out_fn)
    split_num += 1
    print (command)

  return


##
# Main
##
@util.time_dec
def main():
  #print NAME  
  
  # Function calls
  for fn in os.listdir(inp_dir):
    #print fn
    for num in range(700,701):
      #print  'GEN00140{0}'.format(num)
      #print 'GEN00140{0}'.format(num) in fn
      if 'GEN00140{0}'.format(num) in fn and fn[-6:]==".fastq":
        #print fn
        split(os.path.join(inp_dir , fn), fn.replace('.fastq', ''))

  return


if __name__ == '__main__':
  main()
