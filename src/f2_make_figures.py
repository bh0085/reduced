'''
Create figures for reduced library analysis
'''

import os, sys
import pandas as pd
import numpy as np
sys.path.append("/cluster/bh0085")
from mybio import util
from _config import REDUCED_LIB, SEQUENCING_INFO, OUT_PLACE, PRJ_NAME, DATESTR, FIGS_PLACE
import _config

from scipy.stats import entropy
import scipy.stats as stats


import matplotlib.pyplot as plt
#%matplotlib inline

lib_results = pd.read_csv(os.path.join(OUT_PLACE,f"{PRJ_NAME}_aggregate_stats.csv"))

DEFAULT_INP_DIR = _config.OUT_PLACE + 'f1_agg_reports/'
print(DEFAULT_INP_DIR)

if not "__file__" in vars(): __file__="f_test"
NAME = util.get_fn(__file__)
FIG_DIR = os.path.join(FIGS_PLACE, NAME)
util.ensure_dir_exists(FIG_DIR)

agg_results = pd.read_csv(os.path.join(DEFAULT_INP_DIR,f"{PRJ_NAME}_aggregate_stats.csv"))
lib_results = pd.read_csv(os.path.join(DEFAULT_INP_DIR,f"{PRJ_NAME}_lib_stats.csv"))

f,subs = plt.subplots(1,2)
f.set_size_inches(10,4)

plt.sca(subs[0])
ax = plt.gca()
plt.bar(lib_results.index,lib_results.frac_crispr_of_mapped)
ax.set_title("CRISPR cutting efficiency\n(Fraction non-wildtype reads)")
ax.set_ylabel("cutting efficiency")
ax.set_xlabel("library ID")

plt.sca(subs[1])
ax = plt.gca()

plt.bar(lib_results.index, lib_results.deletion_length_r,1, label="all guides")
significant = lib_results.loc[lib_results.deletion_length_pval < .05]
plt.bar(significant.index, significant.deletion_length_r,1,label="significant p<.05")
plt.legend()

#ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis
#ax2.bar(lib_results.index+.5,np.min([-1* np.log10(lib_results.deletion_length_pval),-log10(.001)],.3,color="red")


ax.set_title("CRISPR deletion length spectrum, goodness of fit (pearson R)")
ax.set_ylabel("pearson R")
ax.set_xlabel("library ID")
plt.tight_layout()
plt.savefig("")


f,subs = plt.subplots(1,1)

plt.sca(subs)
x = lib_results.predicted_ins_of_crispr
y = lib_results.frac_ins_of_crispr.fillna(0)


gradient, intercept, r_value, p_value, std_err = sstats.linregress(x,y)
mn=np.min(x)
mx=np.max(x)
x1=np.linspace(mn,mx,500)
y1=gradient*x1+intercept

plt.scatter(x,y,label="guide ID")
plt.plot(x1,y1,'-r',label=f"linear fit\npval = {p_value:0.03}")

plt.savefig(os.path.join(FIG_DIR,"fig1_cutting_efficiency.png"))

ax = plt.gca()
ax.set_xlabel("indelphi predicted insertion fraction")
ax.set_ylabel("observed insertion fraction")
ax.set_title("INSERTION FRACTION\npredicted vs observed for 48 Reduced lib guides")
plt.legend()


plt.savefig(os.path.join(FIG_DIR,"fig2_insertion_fraction.png"))
