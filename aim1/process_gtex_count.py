#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import numpy as np
import pandas as pd
import os
import sys
import pickle
import h5py
import matplotlib.pyplot as plt
import anndata as ad
from scipy.sparse import csr_matrix

'''
Let's assume the GTEx ananlysis is performed using STAR, Ensembl85
And the Exon table is Ensembl 72.

this script is for processing gtex junction counts
'''

# argument
data_folder = '/data/salomonis2/LabFiles/Frank-Li/refactor/data'

# GTEx metadata
sra_table = pd.read_csv(os.path.join(data_folder,'GTEx_SRARunTable.txt'),sep='\t')
sra_table = sra_table.loc[:,['Run','body_site']]   # 24455 non-redundant SRR

# GTEx count
sra_data = pd.read_csv(os.path.join(data_folder,'GTEx_junction_counts.txt'),sep='\t',index_col=0)  # 908631, 30% zero, no missing value
sra_data.columns.name = None
sra_data.index.name = None
sra_data = sra_data.loc[np.logical_not(sra_data.index.duplicated()).tolist(),:]
sra_data.columns = [item.split('_')[0].split(':')[1] for item in sra_data.columns]  # 901299
adata = ad.AnnData(X=csr_matrix(sra_data.values),var=pd.DataFrame(index=sra_data.columns),obs=pd.DataFrame(index=sra_data.index))  

# only need the overlapping samples 
common = list(set(sra_table['Run'].values).intersection(set(adata.var_names)))
sra_table = sra_table.loc[sra_table['Run'].isin(common),:]
srr_to_tissue = sra_table.set_index('Run').squeeze().to_dict()
adata.var['tissue'] = adata.var_names.map(srr_to_tissue).values
print(adata)
adata.write(os.path.join(data_folder,'GTEx_junction_counts.h5ad'))













