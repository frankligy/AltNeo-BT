#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import numpy as np
import pandas as pd
import os
import sys
import pickle
import h5py
import matplotlib.pyplot as plt
import anndata

'''
this script is to query the tumor specificity of the junction
'''

def cruel_tumor_specificity(uid,count):
    adata = anndata.read_h5ad('../data/GTEx_junction_counts.h5ad')
    if uid not in set(adata.obs_names):
        mean_value = 0
    else:
        mean_value = adata.obs.loc[uid,'mean']
    diff = count - mean_value
    if mean_value <= 3 and diff >= 20:
        identity = True
    else:
        identity = False
    return identity













