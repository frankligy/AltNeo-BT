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
from scipy import stats
from scipy.optimize import minimize

'''
calculate neojunction tumor specificity, psi is not informative, count is more informative
'''

data_folder = '/data/salomonis2/LabFiles/Frank-Li/refactor/data'

# # PSI
# with open(os.path.join(data_folder,'gtex_data_df.p'),'rb') as f:
#     sra_data = pickle.load(f)
# sra_data.apply(lambda x:np.ma.masked_invalid(x.values).mean(),axis=1).to_csv('mean_psi.txt',sep='\t')


# # calculate tumor specificity using MLE
# adata = ad.read_h5ad(os.path.join(data_folder,'GTEx_junction_counts.h5ad'))
# adata.obs['mean'] = np.array(adata.X.mean(axis=1)).squeeze()
# adata.obs['std'] = adata.X.toarray().std(axis=1)
# total_count = np.array(adata.X.sum(axis=0)).squeeze() / 1e6
# adata.var['total_count'] = total_count

# def mle_func(parameters,y):
#     sigma = parameters
#     ll = np.sum(stats.halfnorm.logpdf(y,0,sigma))
#     neg_ll = -1 * ll
#     return neg_ll

# query = 'ENSG00000186471:E3.1-E4.1'
# mle_model = minimize(mle_func,np.array([0.2]),args=(adata[query,:].X.toarray().squeeze() / adata.var['total_count'].values,),bounds=((0,1),),method='L-BFGS-B')
# print(mle_model)

# sys.exit('stop')

# calculate tumor specificity using PyMC3
adata = ad.read_h5ad(os.path.join(data_folder,'GTEx_junction_counts.h5ad'))
adata.obs['mean'] = np.array(adata.X.mean(axis=1)).squeeze()
adata.obs['std'] = adata.X.toarray().std(axis=1)
total_count = np.array(adata.X.sum(axis=0)).squeeze() / 1e6
adata.var['total_count'] = total_count

query = 'ENSG00000274565:U0.1_62627849-E1.1_62627009'
y = adata[query,:].X.toarray().squeeze() / adata.var['total_count'].values
print(y)

x = []
for tissue in adata.var['tissue'].unique():
    sub = adata[query,adata.var['tissue']==tissue]
    c = np.count_nonzero(sub.X.toarray())
    x.append(c)
x = np.array(x)

print(x)

import pymc3 as pm
import theano
with pm.Model() as m:
    sigma = pm.Uniform('sigma',lower=0,upper=1)
    nc = pm.HalfNormal('nc',sigma=sigma,observed=y)
    rate = pm.math.sum(nc)/len(y)
    c = pm.Poisson('c',mu=rate,observed=x)

with m:
    step = pm.NUTS()
    trace = pm.sample(500,step=step,return_inferencedata=False,cores=1)

import arviz as az
with m:
    az.plot_trace(trace)
    plt.savefig('arviz.pdf',bbox_inches='tight')
    plt.close()
    df = az.summary(trace,round_to=2)
    print(df)











