'''
Process results of reduced library experiment set
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

# Default params
DEFAULT_INP_DIR = _config.OUT_PLACE + 'e_newgenotype_rev/'

if not "__file__" in vars(): __file__="f_test"
NAME = util.get_fn(__file__)
OUT_DIR = os.path.join(OUT_PLACE, NAME)
util.ensure_dir_exists(OUT_DIR)

results =pd.concat( [
    v for v in SEQUENCING_INFO.Name.apply(
        lambda x:pd.read_csv(
            os.path.join(DEFAULT_INP_DIR, x+"_genotypes_0.csv"),
            index_col=0).assign(sample_name=x)).values]
    ,ignore_index=True)

results["metacat"] = results.apply( lambda x: {"wildtype":"wildtype","del":"del","ins":"ins"}.get(x.Category,"other"),axis=1)
results = results.loc[(results.Length.isna() ) | (results.Length< 15)]
# valid_names = SEQUENCING_INFO.loc[18:77].Name
# results = results.loc[results.sample_name.isin(valid_names)]


# valid_names = SEQUENCING_INFO.Name
# results = results.loc[results.sample_name.isin(valid_names)]
results = results.join(SEQUENCING_INFO.set_index("Name")[["Description"]],on="sample_name")
results = results.join(SEQUENCING_INFO.set_index("Name")[["drug_name","replicate"]],on="sample_name")


#defines a helper function to get groups, return a count of zero if group does not exist
def gg(g, key):
    if key in g.groups: return g.get_group(key)
    else: return pd.DataFrame(pd.Series({"Count":0})).T

    
agg_results = pd.DataFrame()
agg_results["n_ins"] = results.groupby(["sample_name"]).apply(lambda x:gg(x.groupby("Category"),"ins").Count.sum())
agg_results["n_all_del"] = results.groupby(["sample_name"]).apply(lambda x:gg(x.groupby("Category"),"del").Count.sum())
del_results = results.loc[results.Category=="del"]
agg_results["n_mh_del"] = del_results.groupby(["sample_name"]).apply(lambda x:x.loc[x["Microhomology-Based"]=="yes"].Count.sum())
agg_results["n_non_mh_del"] = del_results.groupby(["sample_name"]).apply(lambda x:x.loc[x["Microhomology-Based"]=="no"].Count.sum())

agg_results["n_cas9_total"] = agg_results.n_ins + agg_results.n_all_del
agg_results["n_wildtype"] = results.groupby(["sample_name"]).apply(lambda x:gg(x.groupby("Category"),"wildtype").Count.sum())
agg_results["n_other"] = results.groupby(["sample_name"]).apply(lambda x:gg(x.groupby("metacat"),"other").Count.sum())
agg_results["n_total"] = results.groupby(["sample_name"]).apply(lambda x:x.Count.sum())
agg_results["n_total_mapped"] = results.loc[results.metacat != "other"].groupby("sample_name").Count.sum()


lib_results = pd.DataFrame()
lib_results["n_ins"] = results.groupby(["_Experiment"]).apply(lambda x:gg(x.groupby("Category"),"ins").Count.sum())
lib_results["n_all_del"] = results.groupby(["_Experiment"]).apply(lambda x:gg(x.groupby("Category"),"del").Count.sum())
del_results = results.loc[results.Category=="del"]
lib_results["n_mh_del"] = del_results.groupby(["_Experiment"]).apply(lambda x:x.loc[x["Microhomology-Based"]=="yes"].Count.sum())
lib_results["n_non_mh_del"] = del_results.groupby(["_Experiment"]).apply(lambda x:x.loc[x["Microhomology-Based"]=="no"].Count.sum())

lib_results["n_cas9_total"] = lib_results.n_ins + lib_results.n_all_del
lib_results["n_wildtype"] = results.groupby(["_Experiment"]).apply(lambda x:gg(x.groupby("Category"),"wildtype").Count.sum())
lib_results["n_other"] = results.groupby(["_Experiment"]).apply(lambda x:gg(x.groupby("metacat"),"other").Count.sum())
lib_results["n_total"] = results.groupby(["_Experiment"]).apply(lambda x:x.Count.sum())
lib_results["n_total_mapped"] = results.loc[results.metacat != "other"].groupby("_Experiment").Count.sum()


#annotate results sets with extra information
for results_set in [agg_results,lib_results]:
    results_set["frac_non_mh_del_of_crispr"] = results_set.n_non_mh_del / results_set.n_cas9_total
    results_set["frac_mh_del_of_crispr"] = results_set.n_mh_del / results_set.n_cas9_total
    results_set["frac_del_of_crispr"] = results_set.n_all_del / results_set.n_cas9_total
    results_set["frac_ins_of_crispr"] = results_set.n_ins / results_set.n_cas9_total
    results_set["frac_crispr_of_mapped"] =(results_set.n_all_del +results_set.n_ins) / results_set.n_total_mapped
    
    
    
#load indelphi genotypes for comparison
id_genos = pd.read_csv(os.path.join(OUT_PLACE,"f0_py2_write_indelphi_genotypes/indelphi_genotypes.csv"))

for k,guide in REDUCED_LIB.iterrows():
    for ct in ["mESC"]:
        pred_df = id_genos.loc[(id_genos.celltype==ct) &(id_genos.libid==k) ]
        crispr_results = results.loc[(results._Experiment == k ) & (results.Category.isin(["ins","del"]))]
        fracs = crispr_results.groupby(["Category","Length"]).Count.sum() / crispr_results.Count.sum()
        observed_precision = 1 - entropy(fracs) / np.log(len(fracs))
        lib_results.at[k,"observed_precision"] = observed_precision
        
        overall_precision = 1 - entropy(pred_df['Predicted frequency']) / np.log(len(pred_df))
        highest_fq = max(pred_df['Predicted frequency'])
        highest_del_fq = max(pred_df[pred_df['Category'] == 'del']['Predicted frequency'])/100
        highest_ins_fq = max(pred_df[pred_df['Category'] == 'ins']['Predicted frequency'])/100

        # Outcomes
        ins_fq = sum(pred_df[pred_df['Category'] == 'ins']['Predicted frequency']) / 100
        crit = (pred_df['Category'] == 'del') & (pred_df['Genotype position'] != 'e')
        mhdel_fq = sum(pred_df[crit]['Predicted frequency']) / 100

        crit = (pred_df['Category'] == 'del') & (pred_df['Genotype position'] == 'e')
        nomhdel_fq = sum(pred_df[crit]['Predicted frequency']) / 100

        lib_results.at[k,"predicted_del_of_crispr"] = mhdel_fq + nomhdel_fq
        lib_results.at[k,"predicted_ins_of_crispr"] = ins_fq
        lib_results.at[k,"predicted_mhdel_of_crispr"] = mhdel_fq
        lib_results.at[k,"predicted_non_mhdel_of_crispr"] = nomhdel_fq
        lib_results.at[k,"predicted_precision"] = overall_precision
        
        deletion_lengths = pd.DataFrame(pred_df.loc[pred_df.Category == "del"].groupby("Length")["Predicted frequency"].sum().rename("predicted_freq"))
        for k2,row in deletion_lengths.iterrows():
            deletion_lengths.at[k2,"observed_freq"] = crispr_results.loc[crispr_results.Length==k2].Count.sum()
        deletion_lengths = deletion_lengths.loc[:20].fillna(0)
        
        deletion_length_r,pval = stats.pearsonr(deletion_lengths.predicted_freq,deletion_lengths.observed_freq)
        lib_results.at[k,"deletion_length_r"] = deletion_length_r
        lib_results.at[k,"deletion_length_pval"] = pval


agg_results = agg_results.join(SEQUENCING_INFO.set_index("Name")[["Description","drug_name","replicate"]])
lib_results = lib_results.join(SEQUENCING_INFO.set_index("Name")[["Description","drug_name","replicate"]])

resultsfile = os.path.join(OUT_DIR,f"{PRJ_NAME}_all_results.csv")
results.to_csv(resultsfile)
print(f"wrote all results to {resultsfile}")

#
#note, the transpose of this file gives per-sample results in a useful format for excel users
#we take this transpose when writing to file
#
libfile = os.path.join(OUT_DIR,f"{PRJ_NAME}_library_stats.csv")
lib_results.to_csv(libfile)
print(f"wrote per-library results to {libfile}")

aggfile = os.path.join(OUT_DIR,f"{PRJ_NAME}_aggregate_stats.csv")
agg_results.to_csv(aggfile)
print(f"wrote aggregate results to {aggfile}")






    