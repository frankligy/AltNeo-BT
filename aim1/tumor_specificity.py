#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import numpy as np
import pandas as pd
import os
import sys
import pickle
import h5py
import matplotlib.pyplot as plt

'''
calculate neojunction tumor specificity
'''

data_folder = '/data/salomonis2/LabFiles/Frank-Li/refactor/data'
with open(os.path.join(data_folder,'gtex_data_df.p'),'rb') as f:
    sra_data = pickle.load(f)

sra_data.apply(lambda x:np.ma.masked_invalid(x.values).mean(),axis=1).to_csv('mean_psi.txt',sep='\t')











