#!/Users/ligk2e/opt/anaconda3/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 17 12:27:22 2020

@author: ligk2e
"""
import os
#os.chdir('/Users/ligk2e/Desktop/project_test')
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from yattag import Doc
import h5py
import subprocess

def scratchPlusView1(dataFolder,outFolder):     # directly use pool for multiprocessing

    global dicSRA

    sraTable = pd.read_csv(os.path.join(dataFolder,'GTEx_SRARunTable.txt'),sep='\t')
    sraData = pd.read_csv(os.path.join(dataFolder,'GTEx_EventAnnotation.txt'),sep='\t')
    
    conversion = sraTable[['Run','body_site']]   # this is a complete table
    
    analyzed = sraData.columns.tolist()[11:]
    SRR_ana = [item.split('_')[0] for item in analyzed]    # out data only has 1274 samples
    
    
    
    conversion.index = conversion['Run'].tolist()
    
    new_conversion = conversion.loc[list(set(conversion['Run'].tolist()).intersection(set(SRR_ana)))]  # take intersection
    
    
    
    dicSRA = {}
    for i in range(new_conversion.shape[0]):
        SRR = new_conversion.iloc[i]['Run']
        tissue = new_conversion.iloc[i]['body_site']
        try: 
            dicSRA[tissue].append(SRR)
        except KeyError:
            dicSRA[tissue] = []
            dicSRA[tissue].append(SRR) 
            
    
    
    class Spliter():
        def __init__(self,good):   # good is seperable based on some index
            self.good = good
            self.queue = []  # to store splited data
        
        def split_df(self,n):
            dim  = self.good.shape[0]
            lis = [i for i in range(dim)]
            size = len(lis)//n + 1
            for j in range(0,len(lis),size):
                part_index = lis[j:j+size]
                part_df = self.good.iloc[part_index]
                
                self.queue.append(part_df)

        def split_ndarray(self,n):
            dim = len(self.good)
            lis = [i for i in range(dim)]
            size = len(lis)//n + 1
            for j in range(0,len(lis),size):
                part_ndarray = self.good[j:j+size]
                self.queue.append(part_ndarray)
        
    Spliter1 = Spliter(sraData) 
    Spliter1.split_df(20)  
    liaisonData = Spliter1.queue 
    
    
    import multiprocessing as mp
    p = mp.Pool(processes = 20)
    dicts = p.map(constructDic,liaisonData)
    
    def merge_dicts(dic_list):
    
        result = {}
        for dictionary in dic_list:
            result.update(dictionary)
        return result
    
    dicTissueExp = merge_dicts(dicts)

    with h5py.File(os.path.join(outFolder,'dicTissueExp.hdf5'),'w') as f:
        for event,tissueExp in dicTissueExp.items():
            grp = f.create_group(event)   # /event when grp.name   # i think grp is a group object
            for tissue,expression in tissueExp.items():
                grp.create_dataset(tissue,data=expression)    # /event/tissue  # data has to be a ndarray, not an object

    return dicTissueExp

    # import bz2
    # import _pickle as cpickle
    # import pickle
    # with bz2.BZ2File('./data/dicTissueExp.pbz2','wb') as f1:
    #     cpickle.dump(dicTissueExp,f1)
    # return dicTissueExp

def constructDic(sraData):

    dicTissueExp = {}
    for i in range(sraData.shape[0]):
        print('this is the {0}th run of process{1}'.format(i,os.getpid()))
        event = sraData['UID'].tolist()[i].split('|')[0]
        for tissue,accID in dicSRA.items():
            try: dicTissueExp[event][tissue] = sraData.iloc[i][[accID+'_1.bed' for accID in dicSRA[tissue]]].values.tolist()  # here PSI value will be stored as ndarray
            except KeyError:            
                dicTissueExp[event] = {}
                dicTissueExp[event][tissue] = sraData.iloc[i][[accID+'_1.bed' for accID in dicSRA[tissue]]].values.tolist()
    return dicTissueExp

# def constructHDF5(sraData):

#     with h5py.File('dicTissueExp.hdf5','w') as f:
#         for i in range(sraData.shape[0]):
#             print('this is the {0}th run of process{1}'.format(i,os.getpid()))
#             event = sraData['UID'].tolist()[i].split('|')[0]
#             dic = {}
#             for tissue,accID in dicSRA.items():
#                 try: dic[tissue] = sraData.iloc[i][[accID+'_1.bed' for accID in dicSRA[tissue]]].values  # here PSI value will be stored as ndarray
#                 except KeyError:            
#                     dic[tissue] = {}
#                     dic[tissue] = sraData.iloc[i][[accID+'_1.bed' for accID in dicSRA[tissue]]].values
#                 f.create_dataset(event,data=dic)



    


def loadPlusView():
    from time import process_time
    start = process_time()
    import bz2
    import _pickle as cpickle
    import pickle
    with bz2.BZ2File('./data/dicTissueExp.pbz2','rb') as f1:
         dicTissueExp = cpickle.load(f1)  
    end = process_time()
    
    print('consume {0}'.format(end-start))
    return dicTissueExp



def inspectGTEx(dicTissueExp,event,cutoff,tissue,plot):
    flag=0
    import warnings
    warnings.filterwarnings("ignore")

    if tissue=='all':
        tissueList = []
        summary = []
        try:
            tissueExp = dicTissueExp[event]
        except KeyError:
            summary.append('Don\'t detect expression of {0} in normal tissue'.format(event))
        else:

            for tis,exp in tissueExp.items():
                exp = np.array(exp)   # turn to ndarray, and dtype is object
                exp = exp.astype('float64')
                
                exp[np.isnan(exp)] = 0.0   # because nan just means they don't even have expression

                if exp.size == 0: print('{0} data incomplete'.format(tis))
                elif np.any(exp):   # have non-zero element
                    tissueList.append(tis)
                    hits = sum([True if item > cutoff else False for item in exp]) # how many samples has PSI > 0.1
                    size = exp.size   # how many samples in total

                    out = '{0}:({1}/{2}) has PSI > {3}'.format(tis,hits,size,cutoff)
                    summary.append(out)
                    
                    if plot=='True':
                        fig = plt.figure()
                        plt.bar(np.arange(len(exp)),exp,width=0.2,label=tis)
                        plt.xticks(np.arange(len(exp)),np.arange(len(exp))+1)
                        plt.xlabel('GTEx Samples')
                        plt.ylabel('PSI value')
                        plt.legend()
                        if not os.path.exists(os.path.join(outFolder,'GTEx')): os.makedirs(os.path.join(outFolder,'GTEx'))
                        plt.savefig(os.path.join(outFolder,'GTEx/{0}.svg').format(tis),bbox_inches='tight')
                        plt.close(fig)
                    else: continue
                else: 
                    flag += 1
                    summary.append('No expression in {}'.format(tis))
            summary.append('{0} has no expression in {1} tissue types'.format(event,flag))
            
    else:
        try:
            expression = dicTissueExp[event][tissue]
        except KeyError:
            print('Don\'t detect expression of {0} in normal tissue'.format(event))
        else:
            exp = expression.astype('float64')
            exp[np.isnan(exp)] = 0.0
            if exp.size == 0: print('{0} data incomplete'.format(tissue))
            elif np.any(exp):   # have non-zero element
                
                hits = sum([True if item > cutoff else False for item in exp]) # how many samples has PSI > 0.1
                size = exp.size   # how many samples in total
                print('{0}:({1}/{2}) has PSI > {3}'.format(tissue,hits,size,cutoff))
                
                fig = plt.figure()
                plt.bar(np.arange(len(exp)),exp,width=0.2,label='tissue')
                plt.xticks(np.arange(len(exp)),np.arange(len(exp))+1)
                plt.xlabel('GTEx Samples')
                plt.ylabel('PSI value')
                plt.legend()
                if not os.path.exists('./GTEx'): os.makedirs('./GTEx')
                plt.savefig('./{}.svg'.format(tissue),bbox_inches='tight')
                plt.show()

    return tissueList,summary


def inspectGTExHDF5(event,cutoff,tissue,plot):
    flag=0
    import warnings
    warnings.filterwarnings("ignore")

    if tissue=='all':
        tissueList = []
        summary = []
        with h5py.File(os.path.join(dataFolder,'dicTissueExp.hdf5'),'r') as f:
            try:
                tissueExp = f[event]   # tissueExp is a group object
            except:
                summary.append('Don\'t detect expression of {0} in normal tissue'.format(event))
            else:

                for tis,exp in tissueExp.items():   # group object will have items method
                    exp = np.array(exp)   # turn to ndarray, and dtype is object
                    exp = exp.astype('float64')
                    
                    exp[np.isnan(exp)] = 0.0   # because nan just means they don't even have expression

                    if exp.size == 0: print('{0} data incomplete'.format(tis))
                    elif np.any(exp):   # have non-zero element
                        tissueList.append(tis)
                        hits = sum([True if item > cutoff else False for item in exp]) # how many samples has PSI > 0.1
                        size = exp.size   # how many samples in total

                        out = '{0}:({1}/{2}) has PSI > {3}'.format(tis,hits,size,cutoff)
                        summary.append(out)
                        
                        if plot=='True':
                            fig = plt.figure()
                            plt.bar(np.arange(len(exp)),exp,width=0.2,label=tis)
                            plt.xticks(np.arange(len(exp)),np.arange(len(exp))+1)
                            plt.xlabel('GTEx Samples')
                            plt.ylabel('PSI value')
                            plt.legend()
                            if not os.path.exists(os.path.join(outFolder,'GTEx_{0}'.format(event))): os.makedirs(os.path.join(outFolder,'GTEx_{0}'.format(event)))
                            plt.savefig(os.path.join(outFolder,'GTEx_{0}/{1}.svg').format(event,tis),bbox_inches='tight')
                            plt.close(fig)
                        else: continue
                    else: 
                        flag += 1
                        summary.append('No expression in {}'.format(tis))
                summary.append('{0} has no expression in {1} tissue types'.format(event,flag))
    return tissueList,summary


# html generator

def display_header(event):
    doc, tag, text, line = Doc().ttl()
    with tag('div',klass='header'):
        with tag('div',id='header_title'):
            text('GTEx Viewer Report PSI')
        with tag('div',id='header_event'):
            text('{}'.format(event))
    return doc.getvalue()

def display_summary(tissueList):
    doc, tag, text, line = Doc().ttl()
    with tag('div',klass='summary'):
        with tag('h2'):
            text('Summary')
        with tag('ul'):
            with tag('li'):
                with tag('a',href='#basic1'):
                    text('Basic Summary(Expression)')
            with tag('li'):
                with tag('a',href='#basic2'):
                    text('Basic Summary(No Expression)')

            for i in range(len(tissueList)):
                with tag('li'):
                    with tag('a',href='#M{}'.format(i)):
                        text('{}'.format(tissueList[i]))
                             
    return doc.getvalue()


def display_main(tissueList,summary,cutoff,event):
    doc, tag, text, line = Doc().ttl()
    with tag('div',klass='main'):
        with tag('div',klass='module'):
            with tag('h2',id='basic1'):
                text('Basic Summary(Expression-PSI>{})'.format(cutoff))
            for j in range(len(summary)):
                if not summary[j].startswith('No expression'):
                    with tag('p'):
                        text(summary[j])  
            with tag('h2',id='basic2'):
                text('Basic Summary(No expression)')
            for k in range(len(summary)):
                if summary[k].startswith('No expression'):
                    with tag('p'):
                        text(summary[k]) 
        for i in range(len(tissueList)):
            with tag('div',klass='module'):
                with tag('h2',id='M{}'.format(i)):
                    text('{}'.format(tissueList[i]))
                with tag('p'):
                    doc.stag('img',src='./GTEx_{0}/{1}.svg'.format(event,tissueList[i]))
    return doc.getvalue()

def display_footer():
    doc, tag, text, line = Doc().ttl()
    with tag('div',klass='footer'):
        doc.asis('Author: Guangyuan(Frank) Li <li2g2@mail.uc.edu>')
    return doc.getvalue()

def display_all(tissueList,event,summary,cutoff):
    doc, tag, text, line = Doc().ttl()
    doc.asis('<!DOCTYPE html>')
    with tag('html'):
        with tag('head'):
            line('title','{}_report_PSI'.format(event))
            doc.stag('link',rel='stylesheet',type='text/css',href='style.css')
        with tag('body'):
            doc.asis(display_header(event))
            doc.asis(display_summary(tissueList))
            doc.asis(display_main(tissueList,summary,cutoff,event))
            doc.asis(display_footer())
    return doc.getvalue()

def html_generator(tissueList,event,summary,cutoff):
    with open(os.path.join(outFolder,'test.html'),'w') as f1:
        f1.write(display_all(tissueList,event,summary,cutoff))
  

def usage():


    print('Usage:')
    print('python3 queryGTEx_PSI.py -e KYAT3:ENSG00000137944:E4.1-I4.1_88965413 -c 0.1 -m savage -t all -p True -d /data/salomonis2/LabFiles/Frank-Li/python3/data -o /data/salomonis2/LabFiles/Frank-Li/GTEx_Viewer')
    print('Options:')
    print('-e --event : Splicing event you want to interrogate')
    print('-c --cutoff : cutoff value for PSI')
    print('-m --mode : load off-the-shelf GTEx data or generating it from scratch')
    print('-t --tissue : tissue-specific abundance you want to inspect')
    print('-p --plot : Do you want to plot all tissue-specific expression bar charts?')
    print('-d : data folder path')
    print('-o : output folder path')
    print('-h --help : check help information ')
    print('Author: Guangyuan(Frank) Li <li2g2@mail.uc.edu>, PhD Student, University of Cincinnati, 2020') 

def main(event,mode,cutoff,tissue,plot,dataFolder,outFolder):
    if mode == 'denovo':
        dicTissueExp = scratchPlusView1(dataFolder,outFolder)
        tissueList,summary = inspectGTEx(dicTissueExp,event,cutoff,tissue,plot)
        html_generator(tissueList,event,summary,cutoff)
    if mode == 'savage':
        print('Loding GTEx pre-processed hdf5 file')
        tissueList,summary = inspectGTExHDF5(event,cutoff,tissue,plot)
        html_generator(tissueList,event,summary,cutoff)
    if not os.path.exists(os.path.join(outFolder,'style.css')): subprocess.run(['cp',os.path.join(dataFolder,'style.css'),outFolder])



if __name__ == '__main__':
    import getopt
    import sys
    try:
        options, remainder = getopt.getopt(sys.argv[1:],'he:c:m:t:p:d:o:',['help','event=','cutoff=','mode=','tissue=','plot='])
    except getopt.GetoptError as err:
        print('ERROR:', err)
        usage()
        sys.exit(1)
    for opt, arg in options:
        if opt in ('-e','--event'):
            event = arg
            print('Queried examined splicing event:', arg)
        elif opt in ('-c','--cutoff'):
            #print(arg)
            cutoff = float(arg)
            print('cutoff value is:', arg)
        elif opt in ('-m','--mode'):
            mode = arg
            print('Choose the mode to obtain GTEx data',arg)
        elif opt in ('-t','--tissue'):
            tissue = arg
            print('Tissue I want to inspect:', arg)
        elif opt in ('-p','--plot'):
            plot = arg
            print('Generating plot or not:',arg)
        elif opt in ('-d'):
            dataFolder = arg
            print('data folder is:',arg)   
        elif opt in ('-o'):
            outFolder = arg
            print('output folder is:',arg)         
        elif opt in ('--help','-h'):
            usage() 
            sys.exit()       
    main(event,mode,cutoff,tissue,plot,dataFolder,outFolder)
    
    





















