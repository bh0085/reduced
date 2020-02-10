# 
from _config import SEQUENCING_INFO, REDUCED_LIB, OUT_PLACE


import pandas as pd
import numpy as np
import os, sys, re

import sys, os, datetime, subprocess, imp
sys.path.append('/cluster/bh0085/')
from mybio import util


# Default params
if not "__file__" in vars(): __file__="a_test"
NAME = util.get_fn(__file__)
out_dir = OUT_PLACE + NAME + '/'
util.ensure_dir_exists(out_dir)


##run
# Functions
##
def split(inp_fn, out_nm):
  #print (inp_fn)
  inp_fn_numlines = util.line_count(inp_fn)

  num_splits = 15
  split_size = int(inp_fn_numlines / num_splits)
  if num_splits * split_size < inp_fn_numlines:
    split_size += 1
  while split_size % 4 != 0:
    split_size += 1
  #print (f'Using split size {split_size}')

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
#@util.time_dec
def main(): 
  for k, row in SEQUENCING_INFO.iterrows():
      split(row.fastq_path, row.Name)

  return

if __name__ == '__main__':
  main()
