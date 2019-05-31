# 
from __future__ import division
import _config
import sys, os, fnmatch, datetime, subprocess, imp
sys.path.append('/cluster/mshen/')
import numpy as np
from collections import defaultdict
from mylib import util
import pandas as pd
import matplotlib
matplotlib.use('Pdf')
import matplotlib.pyplot as plt
import seaborn as sns

# Default params
inp_dir = _config.READS_DIR
NAME = util.get_fn(__file__)
out_dir = _config.OUT_PLACE + NAME + '/'
util.ensure_dir_exists(out_dir)

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
  print 'Using split size %s' % (split_size)

  split_num = 0
  for idx in range(1, inp_fn_numlines, split_size):
    start = idx
    end = start + split_size  
    out_fn = out_dir + out_nm + '_%s.fastq' % (split_num)
    command = 'tail -n +%s %s | head -n %s > %s' % (start, inp_fn, end - start, out_fn)
    split_num += 1
    print command

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
