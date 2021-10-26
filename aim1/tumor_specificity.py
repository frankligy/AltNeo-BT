#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import numpy as np
import pandas as pd
import os
import sys
import pickle
import h5py
import matplotlib.pyplot as plt
import anndata as ad
import seaborn as sns

'''
calculate neojunction tumor specificity, psi is not informative, count is more informative
'''

data_folder = '/data/salomonis2/LabFiles/Frank-Li/refactor/data'

# # PSI
# with open(os.path.join(data_folder,'gtex_data_df.p'),'rb') as f:
#     sra_data = pickle.load(f)
# sra_data.apply(lambda x:np.ma.masked_invalid(x.values).mean(),axis=1).to_csv('mean_psi.txt',sep='\t')

# counts
adata = ad.read_h5ad(os.path.join(data_folder,'GTEx_junction_counts.h5ad'))
adata.obs['mean'] = np.array(adata.X.mean(axis=1)).squeeze()
adata.obs['std'] = adata.X.toarray().std(axis=1)
adata.obs.sort_values(by='mean',inplace=True)
adata.obs.to_csv('count_mean_std.txt',sep='\t')
sns.ecdfplot(adata.obs['mean'].values)
plt.savefig('mean_ecdf.pdf',bbox_inches='tight')
plt.close()










