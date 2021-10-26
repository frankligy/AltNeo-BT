'''
Author: Frank Li <li2g2@mail.uc.edu>
July 10th 2020 10:36PM

'''

import pandas as pd
import numpy as np
import os
import ast
import argparse
import ast



def extract(k,HLA,folder,taskName,MHC):
    start = k[0]  # assume it will be 8
    start_filePath = os.path.join(folder,'Neoantigen_{0}_{1}.txt'.format(start,taskName))
    start_data = pd.read_csv(start_filePath,sep='\t')
    start_colName = 'HLAIpresented_dict'
    start_data_base = start_data.drop([start_colName,'HLAIpresented_total','HLAIpresented_count'],axis=1)

    for i in k:
        filePath = os.path.join(folder,'Neoantigen_{0}_{1}.txt'.format(i,taskName))
        data = pd.read_csv(filePath,sep='\t')
        colName = 'HLAIpresented_dict'
        targets = data[colName]   # Series object
        start_data_base = start_data_base.join(pd.Series(targets,name='HLAIpresented_dict_{}mer'.format(i)))

    start_data = combine_all_mers(start_data_base,HLA)
    final_data = count_all_mers(start_data)
    
    final_data.to_csv(os.path.join(folder,'merged_result_{0}.txt'.format(taskName)),sep='\t',index=None)


def combine_all_mers(df,HLA):
    combine_col = []
    for i in range(df.shape[0]):
        combine = []
        col8 = df.iloc[i]['HLAIpresented_dict_8mer']
        col9 = df.iloc[i]['HLAIpresented_dict_9mer']
        col10 = df.iloc[i]['HLAIpresented_dict_10mer']
        col11 = df.iloc[i]['HLAIpresented_dict_11mer']
        for col in [col8,col9,col10,col11]:
            target = col
            if target == 'No candidates': continue
            else:
                target = ast.literal_eval(target)   # a dict
                for hla in HLA:
                    strongBinding = target[hla][0]
                    weakBinding = target[hla][1]
                    combine.extend(strongBinding)
                    combine.extend(weakBinding)
        combine = list(set(combine))
        combine_col.append(combine)
    df['neoantigens'] = combine_col
    return df

def count_all_mers(df):
    count_col = []
    for i in range(df.shape[0]):
        neoantigens = df['neoantigens'].iloc[i]
        count = len(neoantigens)
        count_col.append(count)
    df['neoantigen_counts'] = count_col
    return df



def change_HLA_format(HLA):   #'HLA-A32:01,HLA-A68:01,HLA-B40:01,HLA-B51:01,HLA-C02:02,HLA-C16:01'
    result = HLA.split(',')
    return result



def main(args):

    taskName = args.task
    folder_prefix = args.outdir
    HLA = change_HLA_format(args.HLA)


    folder = os.path.join(folder_prefix, 'resultMHC_{0}'.format(taskName))
    k = [8,9,10,11]
    MHC = 'MHCI'
    extract(k,HLA,folder,taskName,MHC)




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='combine 8mer, 9mer, 10mer, 11mer')
    parser.add_argument('--task',type=str,default=None,help='task name, same as task in mhcPresent.py')
    parser.add_argument('--outdir',type=str,default='.',help='output dicrectory, same as outFolder in mhcPresent.py')
    parser.add_argument('--HLA',type=str,default=None,help='HLA type we want to inspect, same as HLA in mhcPresent.py')
    args=parser.parse_args()
    main(args)

