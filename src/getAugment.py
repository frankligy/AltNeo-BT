'''
Author: Frank Li <li2g2@mail.uc.edu>
July 10th 2020 10:36PM

'''

import pandas as pd
import numpy as np
import os
import ast
import argparse

def extract(k,HLA,folder,taskName,MHC):
    with open(os.path.join(folder,'augment_{0}.fasta'.format(taskName)),'w') as f:
        count = 0
        vault = []
        for i in k:
            filePath = os.path.join(folder,'Neoantigen_{0}_{1}.txt'.format(i,taskName))
            data = pd.read_csv(filePath,sep='\t')
            colName = 'HLAIpresented_dict'
            targets = data[colName]   # Series object
            for target in targets:
                if target == 'No candidates': continue
                else:
                    target = ast.literal_eval(target)   # a dict
                    for hla in HLA:
                        strongBinding = target[hla][0]
                        for sb in strongBinding:
                            if not sb in vault:
                                count += 1
                                f.write('>{0}mer_{1}_strongBinding_{2}\n'.format(i,hla,count))
                                f.write(sb+'\n')
                                vault.append(sb)
                        weakBinding = target[hla][1]
                        for wb in weakBinding:
                            if not wb in vault:
                                count += 1
                                f.write('>{0}mer_{1}_weakBinding_{2}\n'.format(i,hla,count))
                                f.write(wb+'\n')
                                vault.append(wb)


def main(args):
    taskName = args.task
    folder_prefix = args.outdir
    HLA = args.HLA.split(',')
    k = [8,9,10,11]
    folder = os.path.join(folder_prefix, 'resultMHC_{0}'.format(taskName))
    MHC='MHCI'
    extract(k,HLA,folder,taskName,MHC)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='get augmented database for proteomic search')
    parser.add_argument('--task',type=str,default=None,help='task name is the same as taskName in mhcPresent.py')
    parser.add_argument('--outdir',type=str,default=None,help='outdir is the same as outFolder in mhcPresent.py')
    parser.add_argument('--HLA',type=str,default=None,help='HLA you want to inspect, same as -H in mhcPresent.py')
    args=parser.parse_args()
    main(args)

