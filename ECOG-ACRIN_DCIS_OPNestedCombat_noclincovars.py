import pandas as pd
import neuroCombat as nC
from sklearn.preprocessing import LabelEncoder
import matplotlib.pyplot as plt
from scipy.stats import ranksums, ttest_ind, ttest_rel, ks_2samp
import os
import sys
import numpy as np
from itertools import permutations
sys.path.append('/Users/Kalina/Documents/CBIG/Code/opnested-combat-main/Functions')
import OPNestedComBat as nested
from sklearn.cluster import KMeans as kmeans


path = "/Users/Kalina/Documents/CBIG/Project_DCIS_R01/"
codepath = '/Users/Kalina/Documents/CBIG/Code/opnested-combat-main/'
datfile = "first_postcontrast_features_processed_n297_resampv2.xls"
covfile = "first_postcontrast_features_covars_n297_resamp_experiments.xls"

#read in feature data
data_df = pd.read_excel(path+datfile, sheet_name='processed_noHarm_features')
data_df = data_df.rename(columns={"SubjectID": "case"})
data_df.to_csv(codepath+"Results"+"data_df_check_vb.csv")


#read in clinical data
categorical_cols = []
continuous_cols = []

#read in batch effects to control for
batch_df = pd.read_excel(path+covfile, sheet_name='batch_effects')
#batch_df.index = batch_df.case
batch_df = batch_df.drop(labels='case',axis=1)
batch_df.resolution = kmeans(n_clusters=3, random_state=0).fit(batch_df.resolution.to_numpy().reshape(-1,1)).labels_
batch_df["field_strength"] = (batch_df["field_strength"] <= 1.5).astype(int)
batch_list = ['field_strength','resolution']#batch_df.keys().tolist()
covars = batch_df

#remove nans
splitlen = len(batch_list)+len(categorical_cols)+len(continuous_cols)
a = pd.concat([covars,data_df],axis=1).dropna()
a.to_csv(codepath+"Results"+"/a_check_vb.csv")
print(splitlen)
covars = a.iloc[:, :splitlen].reset_index(drop=True)
data_df = a.iloc[:, splitlen :].reset_index(drop=True)
caseno = data_df['case']
data_df = data_df.drop(labels='case',axis=1)
dat = data_df.T.apply(pd.to_numeric)
dat.to_csv(codepath+"Results"+"/dat_final_check_vb.csv")

# # FOR GMM COMBAT VARIANTS:
# # Adding GMM Split to batch effects
gmm_df = nested.GMMSplit(dat, caseno, codepath)
gmm_df.to_csv(codepath+"Results"+"/gmm_check_vb.csv")
#gmm_df_merge = covars_df.merge(gmm_df, right_on='Patient',left_on='resolution')
#gmm_df_merge.to_csv(codepath+"Results"+"/gmm_merge_check_vb.csv")
covars['GMM'] = gmm_df['Grouping'] #gmm_df_merge['Grouping']
covars.to_csv(codepath+"Results"+"/covars_final_check_vb.csv")

print(categorical_cols)
print(continuous_cols)
print(batch_list)

# # EXECUTING OPNESTED+GMM COMBAT
# # Here we add the newly generated GMM grouping to the list of batch variables for harmonization
# batch_list = batch_list + ['GMM']

# EXECUTING OPNESTED-GMM COMBAT
# Here we add the newly generated GMM grouping to the list of categorical variables that will be protected during
# harmonization
categorical_cols = categorical_cols + ['GMM']

# Completing Nested ComBat
output_df = nested.OPNestedComBat(dat, covars, batch_list, codepath+"Results", categorical_cols=categorical_cols,
                                  continuous_cols=continuous_cols)
write_df = pd.concat([caseno, output_df], axis=1) # write results fo file
write_df.to_csv(codepath+'Results'+'/features_DCIS_NestedComBat_resampv2.csv')

#Compute the AD test p-values to measure harmonziation performance
nested.feature_ad(dat.T, output_df, covars, batch_list, codepath+'Results')
# Plot kernel density plots to visualize distributions before and after harmonization
nested.feature_histograms(dat.T, output_df, covars, batch_list, codepath+'Results')
