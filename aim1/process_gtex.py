#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import numpy as np
import pandas as pd
import os
import sys
import pickle
import h5py
import matplotlib.pyplot as plt

'''
Let's assume the GTEx ananlysis is performed using STAR, Ensembl85
And the Exon table is Ensembl 72.
'''

# argument
data_folder = '/data/salomonis2/LabFiles/Frank-Li/refactor/data'

# # GTEx metadata
# sra_table = pd.read_csv(os.path.join(data_folder,'GTEx_SRARunTable.txt'),sep='\t')
# sra_table = sra_table.loc[:,['Run','body_site']]   # 24455 non-redundant SRR

# # GTEx PSI
# sra_data = pd.read_csv(os.path.join(data_folder,'GTEx_EventAnnotation.txt'),sep='\t')
# sra_data = sra_data.iloc[:,11:].set_index(sra_data['UID'].values)  # 14% zero psi, 46% nan
# sra_data.columns = [item.split('_')[0] for item in sra_data.columns]
# sra_data.insert(0,'foreground',[item.split('|')[0] for item in sra_data.index])  # 174139 events

# # only need the overlapping samples
# common = list(set(sra_table['Run'].values).intersection(set(sra_data.columns[1:].values)))
# sra_table = sra_table.loc[sra_table['Run'].isin(common),:]
# # tissue_to_srr = sra_table.groupby('body_site')['Run'].apply(lambda x:x.tolist()).to_dict() # {'adipose:[SRR1,SRR2,..]'}
# srr_to_tissue = sra_table.set_index('Run').squeeze().to_dict()

# # remodel the PSI dataframe
# index_mi = [sra_data.index.values.tolist(),sra_data['foreground'].values.tolist()]  # uid, foreground 
# columns_mi = [sra_data.columns[1:].map(srr_to_tissue).values.tolist(),sra_data.columns[1:]]  # tissue, srr
# data = sra_data.iloc[:,1:].values
# sra_data = pd.DataFrame(data=data,index=pd.MultiIndex.from_arrays(index_mi),columns=pd.MultiIndex.from_arrays(columns_mi))
# with open(os.path.join(data_folder,'gtex_data_df.p'),'wb') as f:
#     pickle.dump(sra_data,f)

# # constrct the dict version, for massive amount of query, it will take 20 mins or so
# with open(os.path.join(data_folder,'gtex_data_df.p'),'rb') as f:
#     sra_data = pickle.load(f)
# tissues = sra_data.columns.levels[0].values
# tmp = sra_data.apply(lambda x:{tissue:x.loc[tissue].values for tissue in tissues},axis=1,result_type='reduce') # index1,index2  --> {tissue:array,}
# tmp.index = tmp.index.get_level_values(level=0)
# dic = tmp.to_dict()
# with open(os.path.join(data_folder,'gtex_data_dict.p'),'wb') as f:
#     pickle.dump(dic,f)

# cut it to each tissue type
for tissue in sra_data.columns.levels[0]:
    sub_sra_data = sra_data.loc[:,tissue] # column is not MultiIndex anymore
    with open(os.path.join(data_folder,'by_tissue','{}.p'.format(tissue)),'wb') as f:
        pickle.dump(sub_sra_data,f)











