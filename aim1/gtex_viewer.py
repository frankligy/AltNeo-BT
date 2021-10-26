#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import numpy as np
import pandas as pd
import os
import sys
import pickle
import h5py
import matplotlib.pyplot as plt

'''
this is gtex viewer
'''

# gtex viewer
def gtex_visual(data_folder,query,out_folder='.'):
    with open(os.path.join(data_folder,'gtex_data_df.p'),'rb') as f:
        sra_data = pickle.load(f)
    if type(query) == int:
        info = sra_data.iloc[query,:]
    else:
        query = np.where(sra_data.index.get_level_values(level=0).values==query)[0][0]
        info = sra_data.iloc[query,:]
    title = info.name[0]
    identifier = query
    n_tissue = len(sra_data.columns.levels[0])
    ncols = 5
    fig,axes = plt.subplots(figsize=(20,20),ncols=ncols,nrows=n_tissue//ncols+1,gridspec_kw={'wspace':0.5,'hspace':0.5})
    for i,ax in enumerate(axes.flat):
        if i < n_tissue:
            tissue = sra_data.columns.levels[0].values[i]
            psi = info.loc[tissue].values
            ax.plot(np.arange(len(psi)),psi,marker='o',markersize=4,color='k',markerfacecolor='r',markeredgecolor='r',linestyle='--')
            ax.set_xticks(np.arange(len(psi)))
            ax.set_xticklabels(['s{}'.format(i) for i in np.arange(len(psi))],fontsize=4,rotation=60)
            ax.set_title(tissue,fontsize=8)
            ax.set_ylim(-0.05,1.05)
            ax.set_ylabel('PSI')
        else:
            ax.axis('off')
            break
    fig.suptitle(title,fontsize=10)
    plt.savefig(os.path.join(out_folder,'gtex_visual_{}.pdf'.format(identifier)),bbox_inches='tight')
    plt.close()

data_folder = '/data/salomonis2/LabFiles/Frank-Li/refactor/data'
gtex_visual(data_folder,'STPG1:ENSG00000001460:E15.2-E27.1|ENSG00000001460:I12.2-E27.1','/data/salomonis2/LabFiles/Frank-Li/refactor/gtex_viewer')











