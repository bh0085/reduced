#!/usr/bin/env python2
import sys,os
sys.path.append('/cluster/mshen/prj/inDelphi')
import inDelphi
from scipy.stats import linregress
import numpy as np
import pandas as pd

sys.path.append("/cluster/bh0085")
from mybio import util
from _config import REDUCED_LIB, OUT_PLACE
import imp

if not "__file__" in vars(): __file__="f_test"
NAME = util.get_fn(__file__)
OUT_DIR = os.path.join(OUT_PLACE, NAME)
util.ensure_dir_exists(OUT_DIR)

all_predictions = pd.DataFrame()
for model in ["mESC", "U2OS"]:
    imp.reload(inDelphi)
    inDelphi.init_model(celltype = model)
    for k,row in REDUCED_LIB.iterrows():
        target_seq = row["Designed sequence (61-bp, cutsite at position 34 by 0-index)"]
        CUTSITE = 34
        pred_df, stats = inDelphi.predict(target_seq,CUTSITE)
        pred_df = inDelphi.add_mhless_genotypes(pred_df,stats)
        pred_df = pred_df.assign(**{"celltype":model,
                                 "libid":k
                                 })
        all_predictions = all_predictions.append(pred_df,ignore_index=True)

all_predictions.to_csv(os.path.join(OUT_DIR,"indelphi_genotypes.csv"),index=False)
