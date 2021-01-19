#!/data/salomonis2/LabFiles/Frank-Li/mhc/test_AltNeo/AltNeo_env/bin/python3.6

import sys
print(sys.path)


"""
Created on Mon Mar 30 17:32:27 2020

@author: ligk2e
"""

import os
import sys
import warnings
warnings.filterwarnings("ignore")
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import bz2
import _pickle as cpickle
import pickle
import re
from time import process_time,sleep
from urllib.error import HTTPError
import requests
import xmltodict
import multiprocessing as mp
import subprocess
from pathos.multiprocessing import ProcessingPool as Pool
import getopt
from decimal import Decimal as D
import shelve
import itertools
import statistics


### IW score module
def PSI(df,cutoff):
    ave_PSI, median_PSI, percentage_PSI = [],[],[]
    print('loading dicTissueExp_psi file')
    with bz2.BZ2File(os.path.join(dataFolder,'dicTissueExp.pbz2'),'rb') as f1:
        dicTissueExp_psi = cpickle.load(f1)
    print('successfully load dicTissueExp_psi file into RAM')
    for i in range(df.shape[0]):
        cumulative = 0   # cumulative PSI values, with the intent of calculating average PSI
        lis = []        # record all PSI value for calcating median
        hits = 0        # how many samples are above threshold
        n = 0             # total amounts of samples
        event = df.iloc[i]['UID'].split('|')[0]
        tissueDic = dicTissueExp_psi[event]
        for tis,exp in tissueDic.items():
            exp = np.array(exp)   # turn to ndarray, and dtype is object
            exp = exp.astype('float64')
            exp[np.isnan(exp)] = 0.0   # because nan just means they don't even have expression
            if tis == 'Cells - Cultured fibroblasts' or tis == 'Cells - Leukemia cell line (CML)' or tis == 'Cells - EBV-transformed lymphocyte':  # these tissue are tumor tissue, should be excluded
                continue
            else: 
                if np.any(exp):   # exist expresson in this tissue type
                    cumulative += sum(exp)   # add all PSI value in this tissue to cumulative
                    lis.extend(exp.tolist())       # inject all samples' PSI to lis
                    hits += sum([True if item>cutoff else False for item in exp])   # calculate how many hits and add it to hits
                    n += exp.size           # calculate how many tissues and add to n
                else:    # all zero
                    lis.extend(exp.tolist())
                    n += exp.size
        # calculate three metrics for each event
        ave = cumulative / n
        median = statistics.median(lis)
        percentage = hits / n 
        # append them to final lists
        ave_PSI.append(ave)
        median_PSI.append(median)
        percentage_PSI.append(percentage)
    return ave_PSI,median_PSI,percentage_PSI

def ReadCounts(df,cutoff):
    ave_counts, median_counts, percentage_counts = [],[],[]
    print('loading dicTissueExp_counts file')
    with bz2.BZ2File(os.path.join(dataFolder,'dicTissueExp_counts.pbz2'),'rb') as f1:
        dicTissueExp_counts = cpickle.load(f1)
    print('successfully load dicTissueExp_counts file into RAM')
    for i in range(df.shape[0]):
        cumulative = 0   # cumulative read counts values, with the intent of calculating average PSI
        lis = []        # record all read counts value for calcating median
        hits = 0        # how many samples are above threshold
        n = 0             # total amounts of samples
        event = df.iloc[i]['UID'].split('|')[0]
        event = ':'.join(event.split(':')[1:])   # trim out the gene symbol, only keep ENSG:E3.4-E5.6
        tissueDic = dicTissueExp_counts[event]
        for tis,exp in tissueDic.items():
            exp = np.array(exp)   # turn to ndarray, and dtype is object
            exp = exp.astype('int')
            if tis == 'Cells - Cultured fibroblasts' or tis == 'Cells - Leukemia cell line (CML)' or tis == 'Cells - EBV-transformed lymphocyte':  # these tissue are tumor tissue, should be excluded
                continue
            else: 
                if np.any(exp):   # exist expresson in this tissue type
                    cumulative += sum(exp)   # add all PSI value in this tissue to cumulative
                    lis.extend(exp.tolist())       # inject all samples' PSI to lis
                    hits += sum([True if item>cutoff else False for item in exp])   # calculate how many hits and add it to hits
                    n += exp.size           # calculate how many tissues and add to n
                else:    # all zero
                    lis.extend(exp.tolist())
                    n += exp.size
        # calculate three metrics for each event
        ave = cumulative / n
        median = statistics.median(lis)
        percentage = hits / n 
        # append them to final lists
        ave_counts.append(ave)
        median_counts.append(median)
        percentage_counts.append(percentage)
    return ave_counts,median_counts,percentage_counts



def scoring_process(df_evaluation,cutoff_PSI,cutoff_readcounts,max_mat_ori,min_mat_ori,max_mat_rescale,min_mat_rescale,leading_eigenvector):
    ave_PSI,median_PSI,percentage_PSI = PSI(df_evaluation,cutoff_PSI)
    ave_counts,median_counts,percentage_counts = ReadCounts(df_evaluation,cutoff_readcounts)
    
    mat_eval = np.column_stack((ave_PSI,median_PSI,percentage_PSI,ave_counts,median_counts,percentage_counts))
    print('shape of mat_eval:',mat_eval.shape)
    print(mat_eval)
    mat_eval_new = np.zeros(mat_eval.shape)

    for j in range(mat_eval.shape[1]):
        for i in range(mat_eval.shape[0]):
            new_ij = core_function(max_mat_ori[j],min_mat_ori[j],max_mat_rescale[j],min_mat_rescale[j],mat_eval[i,j])
            mat_eval_new[i,j] = new_ij
    print('shape of mat_eval_new:',mat_eval_new.shape)
    print(mat_eval_new)
    IWscore = []
    for m in range(df_evaluation.shape[0]):
        score = core_IW(mat_eval_new[m,:],leading_eigenvector)
        inverse_score = (-1.0) * float(score)
        sigmoid_score = sigmoid(inverse_score)
        IWscore.append(sigmoid_score)
    df_evaluation['IWscore'] = IWscore
    return df_evaluation



def wrapper_scoring_process(scoreFile,dataFolder):
    scoring = scoreFile
    # max_mat_ori = [9.98564855e-01,1.00000000e+00,1.00000000e+00,1.24263241e+04,1.02370000e+04,1.00000000e+00]
    # min_mat_ori = [0.00092801,0.0,0.00326531,0.01387755,0.0,0.0]
    # max_mat_rescale = [1.6203541,1.246249,0.98267483,72.70393268,80.23512846,0.77875848]
    # min_mat_rescale = [-1.69933277,-1.43360037,-2.65506607,-0.30237015,-0.30285234,-2.89327914]
    # leading_eigenvector = [0.48264742,0.47174347,0.47383551,0.21692184,0.22945297,0.46934607]
    s = shelve.open(os.path.join(dataFolder,'training_parameters_IW'))
    max_mat_ori = s['max_mat_ori']
    min_mat_ori = s['min_mat_ori']
    max_mat_rescale = s['max_mat_rescale']
    min_mat_rescale = s['min_mat_rescale']
    leading_eigenvector = s['leading_eigenvector']
    s.close()
    df_new = scoring_process(scoring,0.1,3,max_mat_ori,min_mat_ori,max_mat_rescale,min_mat_rescale,leading_eigenvector)
    return df_new

def core_function(max_ori,min_ori,max_re,min_re,old):
    new = (max_re * (old - min_ori) + min_re * (max_ori - old)) / (max_ori - min_ori)
    return new

def core_IW(array,weight):
    IW = np.dot(array,weight)     # vectorization
    return IW

def sigmoid(x):
    x = float(x)
    y = 1/(1+np.exp(-x))
    return y



## main program

class Meta():  #inspect an object: dir(), vars(), instanceName.__dict__, mannually set __repr__ or __str__
    def __init__(self, df):
        
        self.df = df
    
    def __repr__(self):
        return {'df':self.df}
    
    def getPercent(self,dfgroup,dfori,name,write=False):
        dic = {}
        for i in range(dfgroup.shape[0]):
            id_, label = dfgroup['TCGA-ID'].tolist()[i],dfgroup['label'].tolist()[i]
            if label == 'R1-V7':
                try: dic['R1-V7'].append(id_)
                except KeyError:
                    dic['R1-V7'] = []
                    dic['R1-V7'].append(id_)
            elif label == 'Healthy':
                try: dic['Healthy'].append(id_)
                except KeyError:
                    dic['Healthy'] = []
                    dic['Healthy'].append(id_)
        
        num_t = len(dic['R1-V7'])
        num_h = len(dic['Healthy'])
 
        percentArray_h = [] 
        percentArray_t = [] 
        dicPercent = {}        
        for j in range(dfori.shape[0]):
            event = dfori['UID'].tolist()[j]    # next 'UID' was converted to 'UID.1'
            
            allvalue_t = list(map(lambda x:round(x,2),dfori.iloc[j][dic['R1-V7']].tolist()))
            truth_t = [True if item == 0 else False for item in allvalue_t]
            allzero_t = sum(truth_t)   # how many zeros in this cluster
            percent_t = allzero_t/num_t  # percentage of zeros
            percentArray_t.append(percent_t)
                
            #allvalue = dfori[dic['R1-V7']].iloc[j].tolist()
            allvalue_h = list(map(lambda x:round(x,2),dfori.iloc[j][dic['Healthy']].tolist()))
            truth_h = [True if item == 0 else False for item in allvalue_h]
            allzero_h = sum(truth_h)   # how many zeros in this cluster
            percent_h = allzero_h/num_h   # percentage of zeros
            percentArray_h.append(percent_h)
            dicPercent[event]=(percent_t,allvalue_t,percent_h,allvalue_h)
            
        col0,col1,col2,col3 = [],[],[],[]
        for k in range(self.df.shape[0]):
            splice = self.df['UID'].tolist()[k]
            per_t,all_t,per_h,all_h = dicPercent[splice][0],dicPercent[splice][1],dicPercent[splice][2],dicPercent[splice][3]
            col0.append(per_t)
            col1.append(all_t)
            col2.append(per_h)
            col3.append(all_h)
        self.df['tumor_zero_percent'] = col0
        self.df['tumor_distribution'] = col1
        self.df['healthy_zero_percent'] = col2
        self.df['healthy_distribution'] = col3
        if write==True: self.df.to_csv('see{0}.txt'.format(name),sep='\t',index=None)
            
    def retrieveJunctionSite(self,dict_exonCoords,dict_fa):
        '''
        Caveat:
        because of the overlapping issue in exonCoordinates, the resultant junction sequence, its both terminus might have one overhang
        '''
        exam_seq = []
        cond = []
        for i in range(self.df.shape[0]):
            temp = uid(self.df,i)
            EnsID = list(temp.keys())[0].split(':')[1]
            exam_site = list(temp.values())[0][0]
            # back_site = list(temp.values())[0][1]
            exam_site_1 = subexon_tran(exam_site.split('-')[0],EnsID,dict_exonCoords,dict_fa,'site1')
            exam_site_2 = subexon_tran(exam_site.split('-')[1],EnsID,dict_exonCoords,dict_fa,'site2')
            exam_seq_join = ','.join([exam_site_1,exam_site_2])
            if '*' in exam_seq_join: cond.append(False)
            else: cond.append(True)
            exam_seq.append(exam_seq_join)
            # back_site_1 = subexon_tran(back_site.split('-')[0],EnsID,dict_exonCoords,dict_fa,'site1')
            # back_site_2 = subexon_tran(back_site.split('-')[1],EnsID,dict_exonCoords,dict_fa,'site2')
            # back_seq_join = ','.join([back_site_1,back_site_2])
            # back_seq.append(back_seq_join)
            
        self.df['exam_seq'] = exam_seq
        #self.df['back_seq'] = back_seq
        self.df['cond'] = cond
        self.df = self.df[self.df['cond']]
        self.df = self.df.drop(columns=['cond'])
        self.df = self.df.set_index(pd.Index(np.arange(self.df.shape[0])))

        
    def matchWithExonlist(self,df_exonlist,dict_exonCoords):
        col1,col2 = [],[]        
        for i in range(self.df.shape[0]):
            temp=uid(self.df,i)
            EnsID=list(temp.keys())[0].split(':')[1]
            Exons_examined = exon_extract(temp,0,EnsID)
            Exons_back = exon_extract(temp,1,EnsID)
            col1.append(core_match(df_exonlist,dict_exonCoords,EnsID,Exons_examined))
            col2.append(core_match(df_exonlist,dict_exonCoords,EnsID,Exons_back))
        self.df['exam_whole_transcripts'] = col1
        #self.df['back_whole_transcripts'] = col2

    def rescueEvent_1(self):    # second round
        for i in range(self.df.shape[0]):
            temp = uid(self.df,i)
            first = self.df['exam_whole_transcripts'].iloc[i]
            hits = sum([True if item else False for item in first])
            if hits == 0: # first round not able to match, launch second_round
                result = second_round_mhc(temp)
                # print(i)
                # print(result)
                # print(self.df['exam_whole_transcripts'].iloc[i])
                self.df['exam_whole_transcripts'].iloc[i] = list(result)
                

    def rescueEvent_2(self):    # third round
        for i in range(self.df.shape[0]):
            temp = uid(self.df,i)
            firstAndsecond = self.df['exam_whole_transcripts'].iloc[i]
            hits = sum([True if item else False for item in firstAndsecond])
            if hits == 0:  # first and second round still not able to match, launch third_round
                result = third_round_mhc(temp,firstAndsecond)
                self.df['exam_whole_transcripts'].iloc[i] = result

    def getORF(self):
        col = []
        for i in range(self.df.shape[0]):
            whole = self.df['exam_whole_transcripts'].iloc[i]  # ['TTTCGGGAGGCC','TTCCCGGGGGAAA'] OR ['','',''] OR ['intron'] OR ['unrecoverable']
            if whole == ['intron'] or whole == ['unrecoverable']:
                col.append(['None'])   # in ORF column, shown as ['None']
            else:
                hits = sum([True if item else False for item in whole])
                if hits == 0: col.append(['None'])
                else:
                    tempArray = []
                    for transcript in whole:
                        if not transcript: tempArray.append('')
                        else:
                            maxTran = transcript2peptide(transcript)
                            tempArray.append(maxTran)   
                    col.append(tempArray)
        self.df['exam_ORF_tran'] = col                         


        
    def getORFaa(self):

        col = []
        for i in range(self.df.shape[0]):
            ORF = self.df['exam_ORF_tran'].iloc[i]   # ['ATGGTCCCGT','ATGGGGGGGCCAAT'] OR ['None'] OR ['','','']
            if ORF == ['None']: col.append(['None'])
            else:
                tempArray = []
                for transcript in ORF:
                    if not transcript: tempArray.append('')
                    else:
                        maxAA = str(Seq(transcript).translate(to_stop=False))
                        tempArray.append(maxAA)
                col.append(tempArray)
        self.df['exam_ORF_aa'] = col
                    



def third_round_mhc(temp,second):  # after second run
    '''
    In second round
    1. [''] means trans-splicing(non-trailing), novel ordinal, intron retention, they don't have '_' in the exon
    2. ['','','','',...''] means newsplicingsite, trans-splicing(trailing) or Alt5,3 but even trimming the trailing part can not match them to existing ones    
    '''

    EnsGID_this = list(temp.keys())[0].split(':')[1]
    exam1 = list(temp.values())[0][0].split('-')[0]   # E22
    exam2 = list(temp.values())[0][0].split('-')[1]   # ENSG:E31
    if second == [''] and 'ENSG' in exam2:  # trans-splicing(non_trailing)
        EnsGID_trans = exam2.split(':')[0]
        exam2_trans = exam2.split(':')[1]
        full_left = single_left_match(exam1,EnsGID_this)
        full_right = single_right_match(exam2_trans,EnsGID_trans)
        full = cat_left_right_asym(full_left,full_right)
        result = full
    elif 'I' in (exam1 + exam2) and second == ['']:   # intron retention
        result = ['intron']
    elif second == ['']:   # novel ordinal    # E3.4(exam1) - E5.1(exam2)
        full_left = single_left_match(exam1,EnsGID_this)
        full_right = single_right_match(exam2,EnsGID_this)
        full = cat_left_right_sym(full_left,full_right)
        result = full
    else: # ['','',''] 
        result = ['unrecoverable']

    return result

def cat_left_right_sym(full_left,full_right):
    result = []
    for i in range(len(full_left)):
        left = full_left[i]
        right = full_right[i]
        if left and right: result.append(left+right)
        else: result.append('')    # this certain transcript doesn't have left subexon and right exon simutaneously,if don't consider cross-concatenation
    return result        
    
            
def cat_left_right_asym(full_left,full_right):   # will store all possible combinations of former and latter sequence
    result = []
    for i in full_left:
        for j in full_right:
            if i and j: result.append(i + j)   # !!!!!!!! here I don't use ',' for seperating left and right portion
            else: result.append('')
    return result
            
            




def single_left_match(exon,EnsGID):   # give you E22, return all sequence E1,E2,...E22
    global dict_exonCoords
    global dictExonList
    global dict_fa
    
    transcripts = dictExonList[EnsGID]
    result = []
    for tran,item in transcripts.items():
        Exons1 = '|' + exon    #|E22
        Exons2 = exon + '|'    # E22|
        if re.search(rf'{re.escape(Exons1)}',item) or re.search(rf'{re.escape(Exons2)}',item):
            exons = item.split('|')
            dict_judge = {}
            for j in range(len(exons)):
                coords = dict_exonCoords[EnsGID][exons[j]]
                strand = coords[1]
                judge = check_exonlist_general(exons,j,strand)
                dict_judge[exons[j]] = judge  
            bucket = []
            for k in range(len(exons)):
                if not exons[k] == exon:
                    bucket.append(exons[k])
                else:
                    bucket.append(exons[k])
                    break
            full_left = ''
            for m in range(len(bucket)):
                coords = dict_exonCoords[EnsGID][bucket[m]]
                strand = coords[1]
                judge = dict_judge[bucket[m]]
                if strand == '+' and judge:   
                    frag = query_from_dict_fa(dict_fa,coords[2],coords[3],EnsGID,coords[1]) # corresponds to abs_start, abs_end, strand
                elif strand == '+' and not judge:
                    frag = query_from_dict_fa(dict_fa,coords[2],int(coords[3])-1,EnsGID,coords[1]) 
                elif strand == '-' and judge:
                    frag = query_from_dict_fa(dict_fa,coords[2],coords[3],EnsGID,coords[1])
                elif strand == '-' and not judge:
                    frag = query_from_dict_fa(dict_fa,int(coords[2])+1,coords[3],EnsGID,coords[1])  # because of the weird
                    # expression of minus strand, need to draw an illustrator to visulize that.
                full_left += frag
            full_left = full_left.replace('\n','')
            result.append(full_left)
        else:
            result.append('')
            
    return result
        
def single_right_match(exon,EnsGID):   # give you E22, return all sequence E22,E23,.....
    global dict_exonCoords
    global dictExonList
    global dict_fa
    
    transcripts = dictExonList[EnsGID]
    result = []
    for tran,item in transcripts.items():
        Exons1 = '|' + exon    #|E22
        Exons2 = exon + '|'    # E22|
        if re.search(rf'{re.escape(Exons1)}',item) or re.search(rf'{re.escape(Exons2)}',item):
            exons = item.split('|')
            dict_judge = {}
            for j in range(len(exons)):
                coords = dict_exonCoords[EnsGID][exons[j]]
                strand = coords[1]
                judge = check_exonlist_general(exons,j,strand)
                dict_judge[exons[j]] = judge  
            bucket = []
            for k in range(len(exons)):
                if not exons[k] == exon:
                    continue
                else:
                    bucket.extend(exons[k:])
                    break
            full_right = ''
            for m in range(len(bucket)):
                coords = dict_exonCoords[EnsGID][bucket[m]]
                strand = coords[1]
                judge = dict_judge[bucket[m]]
                if strand == '+' and judge:   
                    frag = query_from_dict_fa(dict_fa,coords[2],coords[3],EnsGID,coords[1]) # corresponds to abs_start, abs_end, strand
                elif strand == '+' and not judge:
                    frag = query_from_dict_fa(dict_fa,coords[2],int(coords[3])-1,EnsGID,coords[1]) 
                elif strand == '-' and judge:
                    frag = query_from_dict_fa(dict_fa,coords[2],coords[3],EnsGID,coords[1])
                elif strand == '-' and not judge:
                    frag = query_from_dict_fa(dict_fa,int(coords[2])+1,coords[3],EnsGID,coords[1])  # because of the weird
                    # expression of minus strand, need to draw an illustrator to visulize that.
                full_right += frag
            full_right = full_right.replace('\n','')
            result.append(full_right)
        else:
            result.append('')
            
    return result    

def second_round_mhc(temp):
    EnsID=list(temp.keys())[0].split(':')[1]
    exam1 = list(temp.values())[0][0].split('-')[0]
    exam2 = list(temp.values())[0][0].split('-')[1]    
    if '_' in exam1 and not '_' in exam2:    # mode 1
        exam1_exon = exam1.split('_')[0]
        exam1_coord = exam1.split('_')[1]
        query = exam1_exon + '|' + exam2
        result = second_match(EnsID,query,exam1_coord=exam1_coord)


    if '_' not in exam1 and '_' in exam2:   # mode 2
        exam2_exon = exam2.split('_')[0]
        exam2_coord = exam2.split('_')[1]
        query = exam1 + '|' + exam2_exon
        result = second_match(EnsID,query,exam2_coord=exam2_coord)

        
    if '_' in exam1 and '_' in exam2:
        exam1_exon = exam1.split('_')[0]
        exam1_coord = exam1.split('_')[1]                
        exam2_exon = exam2.split('_')[0]
        exam2_coord = exam2.split('_')[1]
        query = exam1_exon + '|' + exam2_exon
        result = second_match(EnsID,query,exam1_coord=exam1_coord,exam2_coord=exam2_coord)

    if not '_' in exam1 and not '_' in exam2:
        result = ['']         # novel ordinal and intron retention and tran-splicing(non-trailing)
    return result

                  
def second_match(EnsID,query,exam1_coord=False,exam2_coord=False): # dictExonList {EnsGID:{EnsTID:exonlist,EnsTID:exonlist}}
    global dict_exonCoords
    global dictExonList
    global dict_fa
    
    
    
    if exam1_coord==False: mode = 2   # trailing occur in latter one 
    if exam2_coord==False: mode = 1   # trailing occur in former one
    if not exam1_coord==False and not exam2_coord==False: mode = 3    # trailing occur in both ones
    #print(mode)
    exam1 = query.split('|')[0]
    exam2 = query.split('|')[1]
    transcripts = dictExonList[EnsID]
    result = []
    for tran,item in transcripts.items():
        Exons1 = '|' + query
        Exons2 = query + '|'   
        Exons12 = '|'  + query + '|'       
        '''
        This patch is for the a bug like this:
        E8.3|E9.1 will match up to E8.1|E8.2|E8.3|E9.10
        so change to the following:
        1. either match to |E8.3|E9.1|
        2. or you are at the beginning, then you have to match ^E8.3|E9.1|
        3. or you are at the end, then you have to match E8.3|E9.1$
        4. or the transcript is just E8.3|E9.1, then scenario3 still suffice
        '''
        if re.search(rf'{re.escape(Exons12)}',item) or re.search(rf'^{re.escape(Exons2)}',item) or re.search(rf'{re.escape(query)}$',item):
            exons = item.split('|')
            dict_judge = {}
            for j in range(len(exons)):
                coords = dict_exonCoords[EnsID][exons[j]]
                strand = coords[1]
                judge = check_exonlist_general(exons,j,strand)
                dict_judge[exons[j]] = judge
            if mode == 1:
                bucket_left, bucket_right = [],[]
                for i in range(len(exons)):
                    if not exons[i] == exam1:
                        bucket_left.append(exons[i])
                    else: 
                        i += 1
                        bucket_right.extend(exons[i:])
                        break    # till now, we throw the exons before queried one to bucket_left, after queried one to bucket_right
                if bucket_left:
                    full_left = ''
                    for k in range(len(bucket_left)):
                        coords = dict_exonCoords[EnsID][bucket_left[k]]
                        strand = coords[1]
                        judge = dict_judge[bucket_left[k]]
                        if strand == '+' and judge:   
                            frag = query_from_dict_fa(dict_fa,coords[2],coords[3],EnsID,coords[1]) # corresponds to abs_start, abs_end, strand
                        elif strand == '+' and not judge:
                            frag = query_from_dict_fa(dict_fa,coords[2],int(coords[3])-1,EnsID,coords[1]) 
                        elif strand == '-' and judge:
                            frag = query_from_dict_fa(dict_fa,coords[2],coords[3],EnsID,coords[1])
                        elif strand == '-' and not judge:
                            frag = query_from_dict_fa(dict_fa,int(coords[2])+1,coords[3],EnsID,coords[1])  # because of the weird
                            # expression of minus strand, need to draw an illustrator to visulize that.
                        full_left += frag
                    full_left = full_left.replace('\n','')
                else:
                    full_left = ''
                
                
                
                coords_query = dict_exonCoords[EnsID][exam1] 
                strand_query = coords_query[1]
                start = int(coords_query[2])
                judge_query = dict_judge[exam1]
                
                
                if strand_query == '+':
    
                    query_frag = query_from_dict_fa(dict_fa,start,int(exam1_coord),EnsID,strand_query)
                    
                elif strand_query == '-':
                    if not judge_query: start = int(coords_query[2])+1
                    query_frag = query_from_dict_fa(dict_fa,start,int(exam1_coord),EnsID,strand_query)
                

                query_frag = query_frag.replace('\n','')
                
                
                
                if bucket_right:
                    full_right = ''
                    for k in range(len(bucket_right)):
                        coords = dict_exonCoords[EnsID][bucket_right[k]]
                        strand = coords[1]
                        judge = dict_judge[bucket_right[k]]
                        if strand == '+' and judge:   
                            frag = query_from_dict_fa(dict_fa,coords[2],coords[3],EnsID,coords[1]) # corresponds to abs_start, abs_end, strand
                        elif strand == '+' and not judge:
                            frag = query_from_dict_fa(dict_fa,coords[2],int(coords[3])-1,EnsID,coords[1]) 
                        elif strand == '-' and judge:
                            frag = query_from_dict_fa(dict_fa,coords[2],coords[3],EnsID,coords[1])
                        elif strand == '-' and not judge:
                            frag = query_from_dict_fa(dict_fa,int(coords[2])+1,coords[3],EnsID,coords[1])  # because of the weird
                            # expression of minus strand, need to draw an illustrator to visulize that.
                        full_right += frag
                    full_right = full_right.replace('\n','')
                else:
                    full_right = ''
                    
                full = full_left + query_frag + full_right
                result.append(full)
            
            if mode == 2:
                bucket_left, bucket_right = [],[]
                for i in range(len(exons)):
                    if not exons[i] == exam2:
                        bucket_left.append(exons[i])
                    else: 
                        i += 1
                        bucket_right.extend(exons[i:])
                        break    # till now, we throw the exons before queried one to bucket_left, after queried one to bucket_right
                #print(bucket_left,bucket_right)
                if bucket_left:
                    full_left = ''
                    for k in range(len(bucket_left)):
                        coords = dict_exonCoords[EnsID][bucket_left[k]]
                        strand = coords[1]
                        judge = dict_judge[bucket_left[k]]
                        if strand == '+' and judge:   
                            frag = query_from_dict_fa(dict_fa,coords[2],coords[3],EnsID,coords[1]) # corresponds to abs_start, abs_end, strand
                        elif strand == '+' and not judge:
                            frag = query_from_dict_fa(dict_fa,coords[2],int(coords[3])-1,EnsID,coords[1]) 
                        elif strand == '-' and judge:
                            frag = query_from_dict_fa(dict_fa,coords[2],coords[3],EnsID,coords[1])
                        elif strand == '-' and not judge:
                            frag = query_from_dict_fa(dict_fa,int(coords[2])+1,coords[3],EnsID,coords[1])  # because of the weird
                            # expression of minus strand, need to draw an illustrator to visulize that.
                        full_left += frag
                    full_left = full_left.replace('\n','')
                else:
                    full_left = ''
                    
                #print(full_left)
                    
                coords_query = dict_exonCoords[EnsID][exam2] 
                #print(coords_query)
                strand_query = coords_query[1]
                try:
                    judge_query = dict_judge[exam2]
                except KeyError:
                    print(dict_judge,exam2,EnsID,query,tran,file=sys.stderr)
                    raise Exception('check aboved')
                end = int(coords_query[3])
                
                
                if strand_query == '+':
                    if not judge_query: end = int(coords_query[3])-1
                    query_frag = query_from_dict_fa(dict_fa,int(exam2_coord),end,EnsID,strand_query)
                    
                elif strand_query == '-':
                    query_frag = query_from_dict_fa(dict_fa,int(exam2_coord),end,EnsID,strand_query)
                

                query_frag = query_frag.replace('\n','')
                
                if bucket_right:
                    full_right = ''
                    for k in range(len(bucket_right)):
                        coords = dict_exonCoords[EnsID][bucket_right[k]]
                        strand = coords[1]
                        judge = dict_judge[bucket_right[k]]
                        if strand == '+' and judge:   
                            frag = query_from_dict_fa(dict_fa,coords[2],coords[3],EnsID,coords[1]) # corresponds to abs_start, abs_end, strand
                        elif strand == '+' and not judge:
                            frag = query_from_dict_fa(dict_fa,coords[2],int(coords[3])-1,EnsID,coords[1]) 
                        elif strand == '-' and judge:
                            frag = query_from_dict_fa(dict_fa,coords[2],coords[3],EnsID,coords[1])
                        elif strand == '-' and not judge:
                            frag = query_from_dict_fa(dict_fa,int(coords[2])+1,coords[3],EnsID,coords[1])  # because of the weird
                            # expression of minus strand, need to draw an illustrator to visulize that.
                        full_right += frag
                    full_right = full_right.replace('\n','')
                else:
                    full_right = ''
                    
                full = full_left + query_frag + full_right
                #print(full)
                result.append(full)
            
            if mode == 3:
                bucket_left, bucket_right = [],[]
                for i in range(len(exons)):
                    if not exons[i] == exam1:
                        bucket_left.append(exons[i])
                    else: 
                        i += 2
                        bucket_right.extend(exons[i:])
                        break    # till now, we throw the exons before queried one to bucket_left, after queried one to bucket_right
                if bucket_left:
                    full_left = ''
                    for k in range(len(bucket_left)):
                        coords = dict_exonCoords[EnsID][bucket_left[k]]
                        strand = coords[1]
                        judge = dict_judge[bucket_left[k]]
                        if strand == '+' and judge:   
                            frag = query_from_dict_fa(dict_fa,coords[2],coords[3],EnsID,coords[1]) # corresponds to abs_start, abs_end, strand
                        elif strand == '+' and not judge:
                            frag = query_from_dict_fa(dict_fa,coords[2],int(coords[3])-1,EnsID,coords[1]) 
                        elif strand == '-' and judge:
                            frag = query_from_dict_fa(dict_fa,coords[2],coords[3],EnsID,coords[1])
                        elif strand == '-' and not judge:
                            frag = query_from_dict_fa(dict_fa,int(coords[2])+1,coords[3],EnsID,coords[1])  # because of the weird
                            # expression of minus strand, need to draw an illustrator to visulize that.
                        full_left += frag
                    full_left = full_left.replace('\n','')
                else:
                    full_left = ''
                    
                coords_query1 = dict_exonCoords[EnsID][exam1]  
                start1 = int(coords_query1[2])
                strand_query1 = coords_query1[1]
                judge_query1 = dict_judge[exam1]
                if strand_query1 == '+':
                    query_frag1 = query_from_dict_fa(dict_fa,start1,int(exam1_coord),EnsID,strand_query)
                elif strand_query1 == '-':
                    if not judge_query1: start1 = int(coords_query1[2]) + 1
                    query_frag1 = query_from_dict_fa(dict_fa,start1,int(exam1_coord),EnsID,strand_query)
                query_frag1 = query_frag1.replace('\n','')
                    
                coords_query2 = dict_exonCoords[EnsID][exam2] 
                strand_query2 = coords_query[1]
                judge_query2 = dict_judge[exam2]
                end2 = int(coords_query2[3])
                
                if strand_query2 == '+':
                    if not judge_query2: end2 = int(coords_query2[3])-1
                    query_frag2 = query_from_dict_fa(dict_fa,int(exam2_coord),end2,EnsID,strand_query)
                elif strand_query == '-':
                    query_frag2 = query_from_dict_fa(dict_fa,int(exam2_coord),end2,EnsID,strand_query)
                    
                query_frag2 = query_frag2.replace('\n','')
                
                '''
                Remember: the arguments to query_from_dict_fa is very simple, it is just the coord[2] and coord[3] in exonlist file,
                no matter which strand it is on. The positive position of the start and end of a segment.
                
                if judge is false:
                    1. '+': coords[3] - 1 
                    2. '-': coords[2] + 1
                
                
                '''
                
                if bucket_right:
                    full_right = ''
                    for k in range(len(bucket_right)):
                        coords = dict_exonCoords[EnsID][bucket_right[k]]
                        strand = coords[1]
                        judge = dict_judge[bucket_right[k]]
                        if strand == '+' and judge:   
                            frag = query_from_dict_fa(dict_fa,coords[2],coords[3],EnsID,coords[1]) # corresponds to abs_start, abs_end, strand
                        elif strand == '+' and not judge:
                            frag = query_from_dict_fa(dict_fa,coords[2],int(coords[3])-1,EnsID,coords[1]) 
                        elif strand == '-' and judge:
                            frag = query_from_dict_fa(dict_fa,coords[2],coords[3],EnsID,coords[1])
                        elif strand == '-' and not judge:
                            frag = query_from_dict_fa(dict_fa,int(coords[2])+1,coords[3],EnsID,coords[1])  # because of the weird
                            # expression of minus strand, need to draw an illustrator to visulize that.
                        full_right += frag
                    full_right = full_right.replace('\n','')
                else:
                    full_right = ''
                    
                full = full_left + query_frag1 + query_frag2 + full_right
                result.append(full)   
        else:
            result.append('')
    return result                
            
            

class NeoJ(Meta):   
    def __init__(self,df,N):
        super().__init__(df)   # super() means Meta
        self.mer = N    # length of mer that we are aiming to get

    def phaseTranslateJunction(self):
        col0,col1 = [],[]
        for i in range(self.df.shape[0]):
            wholeTrans = list(self.df['exam_whole_transcripts'])[i]   # ['TTTCGGGAGGCC','TTCCCGGGGGAAA'] OR ['','',''] OR ['intron'] OR ['unrecoverable']
            ORFtrans = list(self.df['exam_ORF_tran'])[i]   # ['ATGGTCCCGT','ATGGGGGGGCCAAT'] OR ['None'] OR ['','','']
            junctionOri = list(self.df['exam_seq'])[i]  # have comma to delimit former and latter part
            junction = list(self.df['exam_seq'])[i].replace(',','')

            if wholeTrans == ['intron'] or wholeTrans == ['unrecoverable']: 
                col0.append(['None'])
                col1.append(['None'])
            else: 
                hits_w = sum([True if item else False for item in wholeTrans])
                hits_o = sum([True if item else False for item in ORFtrans])
                if hits_w == 0 or hits_o ==0:
                    col0.append(['None'])
                    col1.append(['None'])
                else:
                    phaseArray,peptideArray = [],[]
                    for j in range(len(wholeTrans)):
                        if not wholeTrans[j]: 
                            phaseArray.append('') 
                            peptideArray.append('') 
                        else:
                            '''
                            To better understand the following process: test this code snippet

                            whole = 'AAAAATTTTTCCCCCGGGGG'
                            orf = 'TTTTTCCCCC'
                            junction = 'CCCC,GGGG'
                            junction_n = junction.replace(',','')

                            startORF = whole.find(orf)
                            endORF = startORF + len(orf)-1
                            startJun = whole.find(junction_n[:-1])
                            endJun = startJun + len(junction_n[:-1]) -1 

                            former = junction.find(',')

                            phase = abs(startJun + former - startORF) % 3 
                            '''
                            startJun = wholeTrans[j].find(junction[:-1]) # trim the right-most overhang for simplicity
                            endJun = startJun + len(junction[:-1]) - 1
                            startORF = wholeTrans[j].find(ORFtrans[j])
                            endORF = startORF + len(ORFtrans[j]) - 1
                            if startJun > endORF or endJun < startORF:
                                phase = '*'  # means junction site will not be translated
                                peptide = '*'
                                phaseArray.append(phase)
                                peptideArray.append(phase)
                                #print('The {}th transcript of {}, even though junction site could match with, but optimal ORF suggests that junction site will not be translated'.format(j,list(self.df['UID'])[i]))
                            else: 
                                former = junctionOri.find(',')  # length of former part, also the index of splicesite
                                phase = abs(startJun + former - startORF) % 3  
                                '''
                                # CCCCC,GGGGG
                                the first triplet: ,GGG, no protruding to the left of comma, phase = 0
                                the first triplet: C,GG, one protruding to the left of comma, phase = 1
                                the first triplet: CC,G, two protruding to the left of comman, phase = 2
                                '''
                                N = self.mer
                                
                                if phase == 0: 
                                    front=0 if former - ((N-1)*3) < 0 else (former - (N-1)*3)
                                    
                                    junctionSlice = junction[front: former + ((N-1)*3)].replace(',','')
                                    #print(junctionSlice)
                                elif phase == 1: 
                                    front=0 if former - ((N-1)*3+1) <0 else (former - ((N-1)*3+1))
                                    junctionSlice = junction[front: former + ((N-1)*3+2)].replace(',','')
                                    
                                elif phase == 2: 
                                    front=0 if former - ((N-1)*3+2) <0 else (former - ((N-1)*3+2))
                                    junctionSlice = junction[front: former + ((N-1)*3+1)].replace(',','')
                        
                    
                                peptide = str(Seq(junctionSlice).translate(to_stop=False))

                                phaseArray.append(phase)
                                peptideArray.append(peptide)
            
                    col0.append(phaseArray) 
                    col1.append(peptideArray)
        self.df['phase'] = col0 # ['',1,'*',2] OR ['None']
        self.df['peptide'] = col1  # ['',STTG,'*',TYCTT] OR ['None']
                
    def getNmer(self):
        col = []
        for i in range(self.df.shape[0]):
            #print(i,'first round-Process{}'.format(os.getpid()))
            peptides = list(self.df['peptide'])[i]   

            if peptides == ['None']: col.append(['MANNUAL'])
            else:
                merArray = []
                for peptide in peptides:
                    if peptide == '': merArray.append('')
                    elif peptide == '*': merArray.append('*')
                    else: 
                        tempArray = extractNmer(peptide,self.mer)   # if no Nmer, it will return a [], otherwise, ['DTDDTY','YTUYDD']
                        merArray.append(tempArray)
                col.append(merArray)  # [ '','*',[], ['MTDJDKFDJF','FJDKFJDF'] ]
        self.df['{0}mer'.format(self.mer)] = col
        
                 
    def mannual(self):
        col = []
        for i in range(self.df.shape[0]):
            #print(i,'second mannual round-process{}'.format(os.getpid()))
            merArray = []
            uid,junction,Nmer = self.df['UID'].tolist()[i],self.df['exam_seq'].tolist()[i],self.df['{0}mer'.format(self.mer)].tolist()[i]
            if Nmer == ['MANNUAL']: 
                EnsGID = uid.split('|')[0].split(':')[1]
                x = uid.split('|')
                try: x[0].split(':')[3]   # these try-excepts aim to consider fusion gene situation
                except IndexError: event = x[0].split(':')[2]
                else: event = str(x[0].split(':')[2])+':'+str(x[0].split(':')[3])     # E22-ENSG:E31  
                try: x[1].split(':')[2]
                except IndexError: backEvent = x[1].split(':')[1]
                else: backEvent = x[1].split(':')[2]   


                if 'ENSG' in event:  # E4.3--ENSG00000267881:E2.1_43874384   # trans-splicing with trailing needs to be rescued
            
                    merArray = tranSplicing(event,EnsGID,junction,self.mer,dictExonList,dict_exonCoords)
                    if merArray == [[]]: merArray = ['transplicing, already checked, query subexon not present in known transcript']

                if re.search(r'I.+_\d+',event) or 'U' in event:  # they belong to NewExon type: including UTR type and blank type(the annotation is blank)
                    # don't infer translation frame for this type, for simplicity
                    junctionIndex = junction.find(',')
                    Nminus1 = self.mer - 1 
                    front = 0 if junctionIndex-(Nminus1*3+2) < 0 else junctionIndex-(Nminus1*3+2)
                    junctionSlice = junction[front:junctionIndex+(Nminus1*3+2)].replace(',','') # no need to worry out of index, slicing operation will automatically change it to 'end' if overflow
                    merArray = dna2aa2mer(junctionSlice,self.mer)   # merArray is a nested list

                if re.search(r'E\d+\.\d+_\d+',event):  # they belong to Alt5, alt3 site that can not be rescued in second round
                
                    merArray = newSplicingSite(event,EnsGID,junction,self.mer,dictExonList,dict_exonCoords)   # nested list
                    if merArray == [[]]: merArray = ['newSplicing, already checked, query subexon not present in known transcript']

                if re.search(r'^I\d+\.\d+-',event) or re.search(r'-I\d+\.\d+$',event): # intron
                    merArray = intron(event,EnsGID,junction,dict_exonCoords,dictExonList,self.mer)  # nested list
                    if merArray == [[]]: merArray = ['intron retention, already checked, either former subexon not present in known transcript or matched transcript is not Ensembl-compatible ID']
                
                if re.search(r'E\d+\.\d+-E\d+\.\d+$',event):   # novel ordinal but not able to be rescued by third round
                    # here we extract the backEvent, which is E24.1-E27.1
                    #print(event)
                    backEvent = backEvent.replace('-','|')  # now it is E24.1|E27.1
                    merArray = novelOrdinal(event,backEvent,EnsGID,junction,dictExonList,dict_exonCoords,self.mer)
                    if merArray ==[[]]: merArray = ['novel ordinal splicing event, already checked, its former part can not match up with any existing transcript either']
                    
            col.append(merArray)
        self.df['mannual'] = col
        

    
    def netMHCresult(self,HLA,pathSoftWare,mode,sb=0.5,wb=2.0):
        col1,col2,col3 = [],[],[]
        for i in range(self.df.shape[0]):
            #print(i)
            merList = self.df['{0}mer'.format(self.mer)].tolist()[i]
            if merList == ['MANNUAL']: merList = self.df['mannual'].tolist()[i]
            #print(merList)
            netMHCpan.pepFile(merList) # get query.pep file
            machine = netMHCpan(os.path.join(outFolder,'resultMHC_{0}'.format(taskName),'temp','query{}.pep'.format(os.getpid())),HLA,pathSoftWare,self.mer,mode,sb,wb)
            dic = machine.seperator()
            if not dic == 'No candidates':
                _, data, length = is_dict_empty(dic)
                data = list(set(data))
                length = len(data)
            else:
                data, length = 'No candidates','No candidates'

            col1.append(dic)
            col2.append(data)
            col3.append(length)
            
        self.df['HLAIpresented_dict'] = col1
        self.df['HLAIpresented_total'] = col2
        self.df['HLAIpresented_count'] = col3

            
    def MHCflurry(self,HLA,mode,sb=50,wb=500):
        HLA = HLA.split(',')
        predictor = Class1PresentationPredictor.load()
        col1,col2,col3 = [],[],[]
        for i in range(self.df.shape[0]):
            merList = self.df['{0}mer'.format(self.mer)].tolist()[i]
            if merList == ['MANNUAL']: merList = self.df['mannual'].tolist()[i]
            merList = pre_process(merList)
            if not merList: 
                col1.append('No candidates')
                col2.append('No candidates')
                col3.append('No candidates')
            else:
                df_result = predictor.predict(peptides=merList,alleles=HLA,verbose=0)
                dic_result = post_process(df_result,HLA,sb,wb)
                cond,data,length = is_dict_empty(dic_result)
                if cond:    
                    col1.append('No candidates')
                    col2.append('No candidates')
                    col3.append('No candidates')
                else: 
                    col1.append(dic_result)
                    col2.append(data)
                    col3.append(length)
            
        self.df['HLAIpresented_dict'] = col1
        self.df['HLAIpresented_total'] = col2
        self.df['HLAIpresented_count'] = col3

    def immunogenecity(self,HLA):
        HLA = HLA.split(',')
        col = []
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        model = dilatedCNN().to(device)
        model.load_state_dict(torch.load('/data/salomonis2/LabFiles/Frank-Li/immunogenecity/transformer/model/dilatedCNN_balance.pth',map_location=torch.device('cpu')))
        for i in range(self.df.shape[0]):
            merList = self.df['HLAIpresented_total'].iloc[i]
            if merList == 'No candidates': col.append('No candidates')
            else:
 
                result = construct_df4deeplearningmodel(merList,HLA,model,device)   # [[merlist1,hla1,diff_value],[merlist1,hla2,diff_value],[]...., [merlist3,hla1,diff_value],.......]

                col.append(result)
        self.df['immunogenecity']=col

def is_dict_empty(dic_result):
    data = []
    for i in dic_result.values():
        for o in i:
            data.extend(o)
    if data: return False, data, len(data)
    else: return True, data, len(data)


def construct_df4deeplearningmodel(merList,HLA,model,device):

    cartesian = list(itertools.product(merList,HLA))   # [(merlist1,hla1),(merlist1,hla2),()...., (merlist3,hla1),.......]
    col1 = [tup[0] for tup in cartesian]
    col2 = [tup[1] for tup in cartesian]
    col3 = [0 for _ in range(len(cartesian))]
    ori = pd.DataFrame({'peptide':col1,'HLA':col2,'immunogenecity':col3})
    scoring_dataset = dataset(ori,hla,dic_inventory)

    scoring_loader = DataLoader(scoring_dataset,batch_size=len(cartesian),shuffle=False,drop_last=True)



    model.eval()
    with torch.no_grad():
        for i,(X,y) in enumerate(scoring_loader):

            x = X.unsqueeze(1).permute(0,1,3,2).to(device)
            y_pred = model(x)
    diff = y_pred[:,1] - y_pred[:,0]
    result = fiddle_result(cartesian,diff)
    return result


def fiddle_result(cartesian,diff):
    diff = diff.numpy()
    result = []
    for i in range(len(cartesian)):
        item = cartesian[i]
        item = list(item)
        item.append(diff[i])
        if diff[i] >=  0:
            result.append(item)
    return result


    




def pre_process(merList):
    result = []
    [result.extend(item) for item in merList if isinstance(item,list)]
    result = list(set(result))
    return result

def post_process(df,hla,sb,wb):
    dic = {}
    for h in hla:
        dic[h] = [[],[]]
    
    for i in range(df.shape[0]):
        peptide = df['peptide'].iloc[i]
        affinity = df['affinity'].iloc[i]
        HLA = df['best_allele'].iloc[i]
        if affinity <= sb:    # strong binder
            dic[HLA][0].append(peptide)
        elif affinity > sb and affinity <= wb:    # weak binder
            dic[HLA][1].append(peptide)
        else:
            continue
        
    return dic

class netMHCpan():
    # ../netMHCpan -p test.pep -BA -xls -a HLA-A01:01,HLA-A02:01 -xlsfile my_NetMHCpan_out.xls
    def __init__(self,intFile,HLA,pathSoftWare,length,mode,sb=0.5,wb=2.0):
        self.intFile=intFile
        self.HLA=HLA
        self.pathSoftWare=pathSoftWare
        self.sb = sb    # float, 0.5
        self.wb = wb    # float, 2.0
        self.length = length
        self.filter = wb  # float, 2.0
        self.mode=mode
    
    @staticmethod
    def pepFile(lis):
        result = []
        [result.extend(item) for item in lis if isinstance(item,list)]
        result = list(set(result))

        with open(os.path.join(outFolder,'resultMHC_{0}'.format(taskName),'temp','query{}.pep'.format(os.getpid())),'w') as f1:
            [f1.write('{}\n'.format(mer)) for mer in result]
        # now we have query.pep file
    
    def seperator(self):
        if self.mode=='MHCI': 
            self.runSoftWareI()
            dic = self.postFileI()
            return dic
        elif self.mode=='MHCII':
            self.runSoftWareII()
            dic = self.postFileII()
            return dic 
            

    
    def runSoftWareI(self):
        with open(os.path.join(outFolder,'resultMHC_{0}'.format(taskName),'temp','resultI{}.txt'.format(os.getpid())),'w') as f3:
            subprocess.run([self.pathSoftWare,'-p',self.intFile, '-BA','-a',self.HLA,'-rth', str(self.sb), '-rlt', str(self.wb), '-l',str(self.length),'-t',str(self.wb)],stdout=f3)        
       # will generate a result file  
    
    def runSoftWareII(self):
        with open(os.path.join(outFolder,'resultMHC_{0}'.format(taskName),'temp','resultII{}.txt'.format(os.getpid())),'w') as f4:
            subprocess.run([self.pathSoftWare,'-f',self.intFile,'-inptype', '1', '-a',self.HLA,'-length',str(self.length)],stdout=f4)
    
    def postFileII(self):
        with open(os.path.join(outFolder,'resultMHC_{0}'.format(taskName),'temp','resultII{}.txt'.format(os.getpid())),'r') as f5,open(os.path.join(outFolder,'resultMHC_{0}'.format(taskName),'temp','resultII_parse{}.txt'.format(os.getpid())),'w') as f6:
            for line in f5:
                if line.startswith('#') or line.startswith('-') or line.strip('\n') == '':continue
                elif re.search(r'^\w+',line): continue
                elif re.search(r'Pos',line): continue
                else: f6.write(line)
        try:df=pd.read_csv(os.path.join(outFolder,'resultMHC_{0}'.format(taskName),'temp','resultII_parse{}.txt'.format(os.getpid())),sep='\s+',header=None,index_col=0,names=[str(i+1) for i in range(11)])
        except pd.errors.EmptyDataError: dic = 'No candidates'
        else:
            hlaAllele = df['2'].tolist()  # HLA Allele
            mer = df['3'].tolist()   # kmer amino acid
            level = df['11'].tolist()  # <=SB or <=WB
            hla = self.HLA
    
            hlaList = hla.split(',')
            hlaNum = len(hlaList) 
            dic = {}   # store all the candidates with binding affinity
            for i in range(hlaNum):
                sb,wb=[],[]   # all strong binding neoantigen, all weak binding neoantigen
                hlaQuery = hlaList[i]
                occurence = [k for k in range(len(hlaAllele)) if hlaAllele[k] == hlaQuery]
                for j in occurence:
                    if level[j]=='<=SB': sb.append(mer[j])
                    elif level[j]=='<=WB':wb.append(mer[j])
                dic[hlaList[i]] = (sb,wb)
            self.neoantigen = dic
        return dic   
    
    
    def postFileI(self):
        # HLA = 'HLA-A01:01,HLA-A03:01,HLA-B07:02,HLA-B27:05,HLA-B58:01'
        with open(os.path.join(outFolder,'resultMHC_{0}'.format(taskName),'temp','resultI{}.txt'.format(os.getpid())),'r') as f1, open(os.path.join(outFolder,'resultMHC_{0}'.format(taskName),'temp','resultI_parse{}.txt'.format(os.getpid())),'w') as f2:
            for line in f1:
                if line.startswith('#') or line.startswith('-') or line.strip('\n') == '': continue
                elif re.search(r'^\w+',line): continue
                elif re.search(r'Pos',line): continue
                else: f2.write(line)
        try:df = pd.read_csv(os.path.join(outFolder,'resultMHC_{0}'.format(taskName),'temp','resultI_parse{}.txt'.format(os.getpid())),sep='\s+', header=None,index_col=0)  
        except pd.errors.EmptyDataError: dic = 'No candidates'   
        else:
            hlaAllele = df[1].tolist()  # HLA Allele
            mer = df[2].tolist()   # kmer amino acid
            level = df[17].tolist()  # SB or WB
            hla = self.HLA
            
            hlaList = hla.split(',')
            hlaNum = len(hlaList) 
            dic = {}   # store all the candidates with binding affinity
            for i in range(hlaNum):
                sb,wb=[],[]   # all strong binding neoantigen, all weak binding neoantigen
                hlaQuery = netMHCpan.hlaType(hlaList[i])  # HLA-A01:01 to HLA-A*01:01
                occurence = [k for k in range(len(hlaAllele)) if hlaAllele[k] == hlaQuery]
                [sb.append(mer[j]) if level[j]=='SB' else wb.append(mer[j]) for j in occurence]
                dic[hlaList[i]] = [sb,wb]
            self.neoantigen = dic
        return dic


            
            
            
            
            
    @staticmethod
    def hlaType(hla):   # convert HLA-A01:01 to HLA-A*01:01
        index1 = re.search(r'^HLA-[A-Z]+',hla).span()[1]   # here, index will be 5, which is '0'
        former,latter = hla[0:index1],hla[index1:]
        hlaNew = former + '*' + latter
        return hlaNew
            

        



    

        
def novelOrdinal(event,backEvent,EnsGID,junction,dictExonList,dict_exonCoords,N): # E25.1-E26.1
    former = event.split('-')[0]  # E25.1
    attrs = dict_exonCoords[EnsGID][former] #  chr, strand, start, end
    strand,start,end = attrs[1], attrs[2], attrs[3]
    breaking = int(end) + 1
    merBucket = []
    allTransDict = dictExonList[EnsGID]
    for tran,exonlist in allTransDict.items():
        if former in exonlist:
            try:tranStartIndex = grabEnsemblTranscriptTable(tran)
            except HTTPError: continue
            else:
                if type(tranStartIndex) == int:
                     if strand == '+':
                         remainder = (int(breaking) - tranStartIndex) % 3   
                         if remainder == 0: 
                             front = 0 if junction.find(',') - ((N-1)*3)<0 else junction.find(',') - ((N-1)*3)
                             junctionSlice = junction[front:junction.find(',')+((N-1)*3)].replace(',','')
                         elif remainder == 1: 
                             front = 0 if junction.find(',') - ((N-1)*3+1)<0 else junction.find(',') - ((N-1)*3+1)
                             junctionSlice = junction[front:junction.find(',')+((N-1)*3+2)].replace(',','')
                         elif remainder == 2: 
                             front=0 if junction.find(',') - ((N-1)*3+2)<0 else junction.find(',') - ((N-1)*3+2)
                             junctionSlice = junction[junction.find(',') - ((N-1)*3+2):junction.find(',')+((N-1)*3+1)].replace(',','')
                     elif strand == '-':

                         remainder = (tranStartIndex - int(breaking)) % 3 
                         if remainder == 0: 
                             front = 0 if junction.find(',') - ((N-1)*3)<0 else junction.find(',') - ((N-1)*3)
                             junctionSlice = junction[front:junction.find(',')+((N-1)*3)].replace(',','')
                         elif remainder == 1: 
                             front = 0 if junction.find(',') - ((N-1)*3+1)<0 else junction.find(',') - ((N-1)*3+1)
                             junctionSlice = junction[front:junction.find(',')+((N-1)*3+2)].replace(',','')
                         elif remainder == 2: 
                             front=0 if junction.find(',') - ((N-1)*3+2)<0 else junction.find(',') - ((N-1)*3+2)
                             junctionSlice = junction[junction.find(',') - ((N-1)*3+2):junction.find(',')+((N-1)*3+1)].replace(',','')
                
                     try:peptide = str(Seq(junctionSlice).translate(to_stop=False))
                     except: print('**** preexisting bugs')
                     else:
                        merArray = extractNmer(peptide,N) 
                        merBucket.append(merArray)
        if merBucket == []: merBucket = [[]]  # means no match for the former for preexisting bugs
           
    return merBucket
                            
def tranSplicing(event,EnsGID,junction,N,dictExonList,dict_exonCoords):  # E4.3_47384738-ENSG00000267881:E2.1 or E4.3-ENSG00000267881:E2.1_43898439
    query = event.split('-')[0]
    formerLength = len(junction.split(',')[0])
    merBucket = []
    try: attrs = dict_exonCoords[EnsGID][query]
    except KeyError: #E6.1_3854692-ENSG00000267881:E2.1
        print('complicated situation(both trailing and transplicing):{0},{1}'.format(EnsGID,event))
        trailing = query.split('_')[1]  #  3854692
        breaking = int(trailing)+1   # the nucleotide immediately after comma
    else:
        strand = attrs[1]
        end = attrs[3]
        breaking = int(end)+1             # the nucleotide immediately after comma
    finally:
        allTranDict = dictExonList[EnsGID]
        for tran,exonlist in allTranDict.items():
            if query in exonlist: 
                try:tranStartIndex = grabEnsemblTranscriptTable(tran)
                except HTTPError: continue   # for instance, BC4389439 it is not a Ensemble-valid transcript ID, just pass this transcript
                else:
                    if type(tranStartIndex) == int:
                        if strand == '+':
                            remainder = (int(breaking) - tranStartIndex) % 3   
                            if remainder == 0: 
                                front = 0 if junction.find(',') - ((N-1)*3)<0 else junction.find(',') - ((N-1)*3)
                                junctionSlice = junction[front:junction.find(',')+((N-1)*3)].replace(',','')
                            elif remainder == 1: 
                                front = 0 if junction.find(',') - ((N-1)*3+1)<0 else junction.find(',') - ((N-1)*3+1)
                                junctionSlice = junction[front:junction.find(',')+((N-1)*3+2)].replace(',','')
                            elif remainder == 2: 
                                front=0 if junction.find(',') - ((N-1)*3+2)<0 else junction.find(',') - ((N-1)*3+2)
                                junctionSlice = junction[junction.find(',') - ((N-1)*3+2):junction.find(',')+((N-1)*3+1)].replace(',','')
                        elif strand == '-':
    
                            remainder = (tranStartIndex - int(breaking)) % 3    # assume tranStartIndex returned by Ensembl is also the position regarding to postive strand,actually doesn't matter, negative is the same for taking remainder
                            if remainder == 0: 
                                front = 0 if junction.find(',') - ((N-1)*3)<0 else junction.find(',') - ((N-1)*3)
                                junctionSlice = junction[front:junction.find(',')+((N-1)*3)].replace(',','')
                            elif remainder == 1: 
                                front = 0 if junction.find(',') - ((N-1)*3+1)<0 else junction.find(',') - ((N-1)*3+1)
                                junctionSlice = junction[front:junction.find(',')+((N-1)*3+2)].replace(',','')
                            elif remainder == 2: 
                                front=0 if junction.find(',') - ((N-1)*3+2)<0 else junction.find(',') - ((N-1)*3+2)
                                junctionSlice = junction[junction.find(',') - ((N-1)*3+2):junction.find(',')+((N-1)*3+1)].replace(',','')
                
                        try:peptide = str(Seq(junctionSlice).translate(to_stop=False))
                        except: print('**** preexistig bugs')
                        else:
                            merArray = extractNmer(peptide,N) 
                            merBucket.append(merArray)
        if merBucket == []: merBucket = [[]]  # means no match for the former subexon or preexisting bugs
    return merBucket

def newSplicingSite(event,EnsGID,junction,N,dictExonList,dict_exonCoords):   # alternative 3' or 5', its naked subexons still doesn't match with any existing one

    try: event.split('-')[0].split('_')[1]   # see if the former one has trailing part E1.2_8483494
    except IndexError: 
        former = event.split('-')[0]  # no trailing part, former part just E1.2, means newSplicingSite is in latter part
        query = former
        try:
            attrs = dict_exonCoords[EnsGID][query]
        except KeyError:
            print('{0} in {1}, preexisting bugs'.format(query,EnsGID))
            raise Exception('Please move out this event:{0} in {1}'.format(query,EnsGID))
        end = attrs[3]
        trailing = int(end) + 1   # the coordinates of latter part position1
    else: 
        query = event.split('-')[0].split('_')[0]  # has trailing part, former part should get rid of trailing part
        trailing = int(event.split('-')[0].split('_')[1]) + 1    # the coordinates of latter part position1
    finally:
        merBucket = []   
        attrs = dict_exonCoords[EnsGID][query] # chr, strand, start, end
        strand = attrs[1]        
        allTransDict = dictExonList[EnsGID] 
        for tran,exonlist in allTransDict.items():
            if query in exonlist: 
                try:tranStartIndex = grabEnsemblTranscriptTable(tran)
                except HTTPError: continue   # for instance, BC4389439 it is not a Ensemble-valid transcript ID, just pass this transcript
                else:
                    if type(tranStartIndex) == int:
                        if strand == '+':
                            remainder = (int(trailing) - tranStartIndex) % 3   
                            if remainder == 0: 
                                front = 0 if junction.find(',') - ((N-1)*3)<0 else junction.find(',') - ((N-1)*3)
                                junctionSlice = junction[front:junction.find(',')+((N-1)*3)].replace(',','')
                            elif remainder == 1: 
                                front = 0 if junction.find(',') - ((N-1)*3+1)<0 else junction.find(',') - ((N-1)*3+1)
                                junctionSlice = junction[front:junction.find(',')+((N-1)*3+2)].replace(',','')
                            elif remainder == 2: 
                                front=0 if junction.find(',') - ((N-1)*3+2)<0 else junction.find(',') - ((N-1)*3+2)
                                junctionSlice = junction[junction.find(',') - ((N-1)*3+2):junction.find(',')+((N-1)*3+1)].replace(',','')
                        elif strand == '-':

                            remainder = (tranStartIndex - int(trailing)) % 3 
                            if remainder == 0: 
                                front = 0 if junction.find(',') - ((N-1)*3)<0 else junction.find(',') - ((N-1)*3)
                                junctionSlice = junction[front:junction.find(',')+((N-1)*3)].replace(',','')
                            elif remainder == 1: 
                                front = 0 if junction.find(',') - ((N-1)*3+1)<0 else junction.find(',') - ((N-1)*3+1)
                                junctionSlice = junction[front:junction.find(',')+((N-1)*3+2)].replace(',','')
                            elif remainder == 2: 
                                front=0 if junction.find(',') - ((N-1)*3+2)<0 else junction.find(',') - ((N-1)*3+2)
                                junctionSlice = junction[junction.find(',') - ((N-1)*3+2):junction.find(',')+((N-1)*3+1)].replace(',','')
                
                        try:peptide = str(Seq(junctionSlice).translate(to_stop=False))
                        except:print('**** preexising bugs')
                        else:
                            merArray = extractNmer(peptide,N) 
                            merBucket.append(merArray)
        if merBucket == []: merBucket = [[]]  # means no match for the former subexon.
           
    return merBucket



                    
def dna2aa2mer(dna,N):
    manner1,manner2,manner3 = dna[0:],dna[1:],dna[2:]
    merBucket = []
    for manner in [manner1,manner2,manner3]:
        #print(manner)
        try:
            aa = str(Seq(manner).translate(to_stop=False))
        except: 
            print('There are *** in the junction site, previous bugs')
        else:    
            mer = extractNmer(aa,N)
            merBucket.append(mer)
    if merBucket == []: merBucket = [[]]
    return merBucket # a nested list of 9mer
              
                
             
def extractNmer(peptide,N):  # it already considers '' and '*'
    starIndex = peptide.find('*') 
   # print(starIndex)
    merArray = []
    if starIndex == -1 and len(peptide) >= N:
        for i in range(0,len(peptide)-N+1,1):
            mer = peptide[i:i+N]
            merArray.append(mer)
    if starIndex >= N:
        peptideTrun = peptide[:starIndex]
        #print(peptideTrun)
        for j in range(0,len(peptideTrun)-N+1,1):
            mer = peptideTrun[j:j+N]
            merArray.append(mer)
    return merArray
                
        
        

def transcript2peptide(cdna_sequence):   # actually to ORF
    reading_manners = []
    reading_manners.append(cdna_sequence[0:])
    reading_manners.append(cdna_sequence[1:])
    reading_manners.append(cdna_sequence[2:])
    frag_comp_array = []
    for manner in reading_manners:       
        pos = []
        for m in re.finditer(r'(TAA|TGA|TAG)',manner):   # for multiple instances
            if m.start() % 3 == 0:
                pos.append(m.start())
        if pos == []:
            frag_comp_array.extend(rescue_position(pos,manner))
        else:
            frag_array,last_seq = pos_to_frags(pos,manner)
            for frag in frag_array:
                if 'ATG' not in frag or len(frag) == 0:
                    continue
                else:
                    for n in re.finditer(r'ATG',frag):
                        if (len(frag) - n.start()) % 3 == 0:
                            frag_comp = frag[n.start():]
                            frag_comp_array.append(frag_comp)
                            break   # might have multiple 'ATG' so it is necessary to break when find first 'ATG'
                        else:
                            continue
        # process last_seq:
            for n in re.finditer(r'ATG',last_seq):
                if n.start() % 3 == 0:
                    last_frag = last_seq[n.start():]
                    protruding = len(last_frag) % 3
                    end = -1 - protruding + 1   # python end exclusive, so + 1
                    last_frag_real = last_frag[:end]
                    frag_comp_array.append(last_frag_real)
    #######################  # We think if you only has longer length(0-7) but add_score is not higher than original one, you are FAlSE
    max_seq = ''
    max_length = 0
    max_item_score = 0
    for item in frag_comp_array:
        temp1 = len(item)
        if temp1==0: continue
        else:
            add_score = score_GC(item) + score_coding_bias(item)
            if (temp1 - max_length) >= 8:
                max_length = temp1
                max_item_score = add_score
                max_seq = item
            elif (temp1 - max_length) >= 0 and (temp1 - max_length) < 8:
                if add_score >= max_item_score:
                    max_length = temp1
                    max_item_score = add_score
                    max_seq = item
#           else:
#                print('equal length but less likely to be a true ORF or longer length but less likely to be a true ORF',add_score,max_item_score) 
    max_seq_tran = max_seq
    return max_seq_tran

def rescue_position(pos,manner):
    for m in re.finditer(r'ATG',manner):
        if m.start() % 3 ==0:
            span = len(manner) - m.start()
            protruding = span % 3
            end = -1 - protruding + 1
            frag = manner[m.start():end]
            pos.append(frag)
    return pos   # pos is actually a list of peptides
            


def pos_to_frags(pos,sequence):
    frag_array = []
    if pos:       
        frag_array.append(sequence[0:pos[0]])
        i = 0
        while i < len(pos)-1:
            frag_array.append(sequence[pos[i]+3:pos[i+1]])
            i += 1
        last_seq = sequence[pos[-1]+3:]
#       I think following original if condition is not correct, we should keep the last_sequence
#        if not any(codon in last_seq for codon in ['TAA','TAG','TGA']):
#            frag_array.append(sequence[pos[-1]+3:])
    return frag_array, last_seq  # last_seq need special care


def score_GC(sequence):
    GC_content = 0
    length_seq = len(sequence)
    for nt in sequence:
        if nt == 'G' or nt == 'C':
            GC_content += 1
    try:
        GC_percent = GC_content / length_seq
    except:
        print('here it the seq:',sequence)
        raise Exception
    return GC_percent
            
def score_coding_bias(sequence):
    # coding frequency table is from GenScript webpage
    usage_dict = {'TTT':16.9,'TTC':20.4,'TTA':7.2,'TTG':12.6,'TAT':12.0,'TAC':15.6,'TAA':0.7,'TAG':0.5,
                  'CTT':12.8,'CTC':19.4,'CTA':6.9,'CTG':40.3,'CAT':10.4,'CAC':14.9,'CAA':11.8,'CAG':34.6,
                  'ATT':15.7,'ATC':21.4,'ATA':7.1,'ATG':22.3,'AAT':16.7,'AAC':19.5,'AAA':24.0,'AAG':32.9,
                  'GTT':10.9,'GTC':14.6,'GTA':7.0,'GTG':28.9,'GAT':22.3,'GAC':26.0,'GAA':29.0,'GAG':40.8,
                  'TCT':14.6,'TCC':17.4,'TCA':11.7,'TCG':4.5,'TGT':9.9,'TGC':12.2,'TGA':1.3,'TGG':12.8,
                  'CCT':17.3,'CCC':20.0,'CCA':16.7,'CCG':7.0,'CGT':4.7,'CGC':10.9,'CGA':6.3,'CGG':11.9,
                  'ACT':12.8,'ACC':19.2,'ACA':14.8,'ACG':6.2,'AGT':11.9,'AGC':19.4,'AGA':11.5,'AGG':11.4,
                  'GCT':18.6,'GCC':28.5,'GCA':16.0,'GCG':7.6,'GGT':10.8,'GGC':22.8,'GGA':16.3,'GGG':16.4} 
    # do a normaliztion for each triplet, then for all the triplet's sum, divided by the number of triplet
    min_freq = 4.5
    max_freq = 40.8
    norm_usage_dict = {}
    for codon,freq in usage_dict.items():
        norm_usage_dict[codon] = float((D(freq) - D(min_freq)) / (D(max_freq) - D(min_freq)))        
    length_seq = len(sequence)
    num_triplet = length_seq/3
    i = 0   
    score = 0
    while i < length_seq - 2:
        triplet = sequence[i:i+3:1]
        score_tri = norm_usage_dict[triplet]
        score += score_tri
        i += 3
    score_bias = score/num_triplet # scale by the number of triplet in the sequence
    return score_bias  

def exon_extract(temp,pos,EnsID):
    Exons = list(temp.values())[0][pos].split('-')[0] + '|' + list(temp.values())[0][pos].split('-')[1]
    return Exons

def core_match(df_exonlist,dict_exonCoords,EnsID,Exons):
   
    try:
        df_certain = df_exonlist[df_exonlist['EnsGID'] == EnsID]
    except: full_transcript_store = []  # EnsGID is absent in df_exonlist
    full_transcript_store = []
    for item in list(df_certain['Exons']):
        full_transcript=''
        if Exons in item:
            Exonlist = item.split('|')
            for j in range(len(Exonlist)):
                coords = dict_exonCoords[EnsID][Exonlist[j]]
                strand = coords[1]
                judge = check_exonlist_general(Exonlist,j,strand)
                if strand == '+' and judge:   
                    frag = query_from_dict_fa(dict_fa,coords[2],coords[3],EnsID,coords[1]) # corresponds to abs_start, abs_end, strand
                elif strand == '+' and not judge:
                    frag = query_from_dict_fa(dict_fa,coords[2],int(coords[3])-1,EnsID,coords[1]) 
                elif strand == '-' and judge:
                    frag = query_from_dict_fa(dict_fa,coords[2],coords[3],EnsID,coords[1])
                elif strand == '-' and not judge:
                    frag = query_from_dict_fa(dict_fa,int(coords[2])+1,coords[3],EnsID,coords[1])  # because of the weird
                    # expression of minus strand, need to draw an illustration to visulize that.
                full_transcript += frag
            full_transcript = full_transcript.replace('\n','')
            full_transcript_store.append(full_transcript)   
        else: 
            full_transcript_store.append('')
    return full_transcript_store  # ['','ATTTTT','TTTGGCC'], # [] if EnsGID is not present in exonlist

def check_exonlist_general(exonlist,index,strand):
    dict = {}
    for subexon in exonlist:
        exon_num = subexon.split('.')[0]
        subexon_num = subexon.split('.')[1]
        if exon_num in dict:
            dict[exon_num].append(subexon_num)
        else:
            dict[exon_num] = []
            dict[exon_num].append(subexon_num)  # E14 > 1,2,4,5
    # check
    query = exonlist[index]
    query_exon_num = query.split('.')[0]   #E14.1
    query_subexon_num = int(query.split('.')[1])   #2, it is a int
    if strand == '+':
        if str(query_subexon_num + 1) in dict[query_exon_num]: return False
        else: return True
    else:
        if str(query_subexon_num + 1) in dict[query_exon_num]: return False
        else: return True
        

        
def neoJunction_testMethod_noGTEx(df):
    return df        
        
def check_GTEx(df,cutoff_PSI,cutoff_sampleRatio,cutoff_tissueRatio):
    col = []
    for i in range(df.shape[0]):
        UID = df.iloc[i]['UID']
        event = UID.split('|')[0]       # foreground event
        try:
            tissueExp = dicTissueExp[event]  # {heart:[],brain:[]}   # values are ndarray
        except KeyError:
            cond = True
            col.append(cond)
        else:
            tissueCounter = 0
            for tis,exp in tissueExp.items():
                if tis == 'Cells - Cultured fibroblasts' or tis == 'Cells - Leukemia cell line (CML)' or tis == 'Cells - EBV-transformed lymphocyte':  # these tissue are tumor tissue, should be excluded
                    continue
                else:
    
                    exp = exp.astype('float64')
                    exp[np.isnan(exp)] = 0.0   # nan means can not detect the gene expression
                    hits = sum([True if i > cutoff_PSI else False for i in exp])   # in a tissue, how many samples have PSI > cutoff value
                    total = exp.size    # how many samples for each tissue type
                    sampleRatio = hits/total    # percentage of sampels that are expressing this event
                    if sampleRatio > cutoff_sampleRatio: tissueCounter += 1   # this tissue is expressiing this event
            tissueRatio = tissueCounter/51    # 51 tissue types in total,excluding three cancer cell lines
            if tissueRatio > cutoff_tissueRatio:
                cond = False
                col.append(cond)
            else: 
                cond = True
                col.append(cond)
    df['cond'] = col
    new_df = df[df['cond']]
    new_df = new_df.drop(columns = ['cond'])
    new_df = new_df.set_index(pd.Index(np.arange(new_df.shape[0])))
    return new_df   

def check_GTEx_PSI_pbz2(df,cutoff_PSI,cutoff_sampleRatio,cutoff_tissueRatio,taskName,outFolder,dataFolder):
    with bz2.BZ2File(os.path.join(dataFolder,'dicTissueExp.pbz2'),'rb') as f1:
        dicTissueExp = cpickle.load(f1)  
    col = []
    for i in range(df.shape[0]):
        #print(i)
        UID = df.iloc[i]['UID']
        event = UID.split('|')[0]       # 0 means foreground event, 1 means background event
        try:
            tissueExp = dicTissueExp[event]  # {heart:[],brain:[]}   # values are ndarray
        except KeyError:
            cond = 'Unknown'
            col.append(cond)
        else:
            tissueCounter = 0
            for tis,exp in tissueExp.items():

                if tis == 'Cells - Cultured fibroblasts' or tis == 'Cells - Leukemia cell line (CML)' or tis == 'Cells - EBV-transformed lymphocyte':  # these tissue are tumor tissue, should be excluded
                    continue
                else:
                    exp = exp.astype('float64')
                    exp[np.isnan(exp)] = 0.0   # nan means can not detect the gene expression
                    hits = sum([True if i > cutoff_PSI else False for i in exp])   # in a tissue, how many samples have PSI > cutoff value
                    total = exp.size    # how many samples for each tissue type
                    sampleRatio = hits/total    # percentage of sampels that are expressing this event
                    if sampleRatio > cutoff_sampleRatio: tissueCounter += 1   # this tissue is expressiing this event
            tissueRatio = tissueCounter/51    # 51 tissue types in total, excluding 3 cancer cell lines
            if tissueRatio > cutoff_tissueRatio:
                cond = 'False'
                col.append(cond)
            else: 
                cond = 'True'
                col.append(cond)
    df['check'] = col


    df_train = df.loc[df['check']!='Unknown']
    df_train = df_train.set_index(pd.Index(np.arange(df_train.shape[0])))
    df_train = wrapper_scoring_process(df_train,dataFolder)
    df_train.to_csv(os.path.join(outFolder,'resultMHC_{0}'.format(taskName),'df_train_{0}.txt'.format(taskName)),sep='\t',index=None)  # true and false = train

    df_final = df.loc[df['check']!='False']
    df_final = df_final.set_index(pd.Index(np.arange(df_final.shape[0])))
    df_final.to_csv(os.path.join(outFolder,'resultMHC_{0}'.format(taskName),'df_final_{0}.txt'.format(taskName)),sep='\t',index=None)   # true and unknown = final

    df_true = df.loc[df['check'] == 'True']
    df_true = df_true.set_index(pd.Index(np.arange(df_true.shape[0])))
    df_true.to_csv(os.path.join(outFolder,'resultMHC_{0}'.format(taskName),'df_true_{0}.txt'.format(taskName)),sep='\t',index=None)   # true

    df_unknown = df.loc[df['check'] == 'Unknown']
    df_unknown = df_unknown.set_index(pd.Index(np.arange(df_unknown.shape[0])))
    df_unknown.to_csv(os.path.join(outFolder,'resultMHC_{0}'.format(taskName),'df_unknown_{0}.txt'.format(taskName)),sep='\t',index=None)  # unknown

    df_false = df.loc[df['check'] == 'False']
    df_false = df_false.set_index(pd.Index(np.arange(df_false.shape[0])))
    df_false.to_csv(os.path.join(outFolder,'resultMHC_{0}'.format(taskName),'df_false_{0}.txt'.format(taskName)),sep='\t',index=None)     # false

    df.to_csv(os.path.join(outFolder,'resultMHC_{0}'.format(taskName),'df_whole_{0}.txt'.format(taskName)),sep='\t',index=None)  # true + false + unknown
    return df_final    
        
    
def uid(df, i):
    uid = list(df['UID'])[i]       
    gene = uid.split(':')[0]
    dict = {}
    gene = gene + ':' + uid.split('|')[1].split(':')[0] # slicing the ENSG in background event
    x = uid.split('|')
    try: x[0].split(':')[3]
    except IndexError: event = x[0].split(':')[2]
    else: event = str(x[0].split(':')[2])+':'+str(x[0].split(':')[3])
    finally: dict[gene] = [event]
    
    try: x[1].split(':')[2]
    except IndexError: backEvent = x[1].split(':')[1]
    else: backEvent = str(x[1].split(':')[1])+':'+str(x[1].split(':')[2])
    finally: dict[gene].append(backEvent)

    #{'gene:ENSid':[E22-E33,E34-E56]}
    # if fusion gene: E22-ENSG:E31
    return dict

def subexon_tran(subexon,EnsID,dict_exonCoords,dict_fa,flag): # flag means if it is site_1 or site_2
    try:
        attrs = dict_exonCoords[EnsID][subexon]   #  !!!, this exclamation symbol is for following code to track
        exon_seq = query_from_dict_fa(dict_fa,attrs[2],attrs[3],EnsID,attrs[1])  
    except KeyError:
        if ':' in subexon:   #fusion gene
            fusionGeneEnsID = subexon.split(':')[0] # this kind of subexon must be site2
            fusionGeneExon = subexon.split(':')[1]
            if  '_' in fusionGeneExon:   # ENSG:E2.1_473843893894
                suffix = fusionGeneExon.split('_')[1]
                subexon = fusionGeneExon.split('_')[0]
                attrs = dict_exonCoords[fusionGeneEnsID][subexon]
                exon_seq = query_from_dict_fa(dict_fa,suffix,attrs[3],fusionGeneEnsID,attrs[1])
            else:    # ENSG:E2.1
                try:
                    attrs = dict_exonCoords[fusionGeneEnsID][fusionGeneExon]
                except KeyError:
                    exon_seq = '***********************'
                    print('{0} does not include in {1} exonlists'.format(fusionGeneExon,fusionGeneEnsID))
                else:   
                    exon_seq = query_from_dict_fa(dict_fa,attrs[2],attrs[3],fusionGeneEnsID,attrs[1])
        else:
            try:   #E2.1_67878789798
                suffix = subexon.split('_')[1]
            except IndexError:  # it means it throw a keyerror in last exclamation symbol, means E2.1 doesn't include in exonlist because it doesn't have trailing part either
                exon_seq = '***********************'   
                print('{0} does not include in {1} exonlists'.format(subexon,EnsID))
            else:
                subexon = subexon.split('_')[0]
                try:
                    attrs = dict_exonCoords[EnsID][subexon]
                except KeyError:
                    #print('{0} observes an UTR event {1}'.format(EnsID,subexon))
                    chrUTR,strandUTR = utrAttrs(EnsID,dict_exonCoords)

                    exon_seq = utrJunction(suffix,EnsID,strandUTR,chrUTR,flag)  

                else:
                    if flag == 'site2':           
                        exon_seq = query_from_dict_fa(dict_fa,suffix,attrs[3],EnsID,attrs[1])  # chr,strand, start,end
                    elif flag == 'site1':
                        exon_seq = query_from_dict_fa(dict_fa,attrs[2],suffix,EnsID,attrs[1])
    return exon_seq

def utrJunction(site,EnsGID,strand,chr_,flag):  # U0.1_438493849, here 438493849 means the site
    if flag == 'site1' and strand == '+':  # U0.1_438493849 - E2.2
        otherSite = int(site) - 100 + 1   # extract UTR with length = 100
        exon_seq = retrieveSeqFromUCSCapi(chr_,int(otherSite),int(site))
    elif flag == 'site1' and strand == '-':    # 438493849 is the coordinates in forward strand
        otherSite = int(site) + 100 - 1 
        #exon_seq = query_from_dict_fa(dict_fa,site,otherSite,EnsGID,strand)    # site, otherSite must be coordinates in forward strand, strand argument will handle the conversion automatically
        exon_seq = retrieveSeqFromUCSCapi(chr_,int(site),int(otherSite))
        exon_seq = str(Seq(exon_seq).reverse_complement())
    elif flag == 'site2' and strand == '+':  # E5.3 - U5.4_48374838
        otherSite = int(site) + 100 -1
        exon_seq = retrieveSeqFromUCSCapi(chr_,int(site),int(otherSite))
    elif flag == 'site2' and strand == '-':
        otherSite = int(site) - 100 + 1
        #print(EnsGID,chr_,site,otherSite)
        exon_seq = retrieveSeqFromUCSCapi(chr_,int(otherSite),int(site))
        exon_seq = str(Seq(exon_seq).reverse_complement())
    return exon_seq

def utrAttrs(EnsID,dict_exonCoords):  # try to get U0.1's attribute, but dict_exonCoords doesn't have, so we just wanna get the first entry for its EnsGID
    exonDict = dict_exonCoords[EnsID] 
    attrs = next(iter(exonDict.values()))
    chr_,strand = attrs[0],attrs[1]
    return chr_,strand



def retrieveSeqFromUCSCapi(chr_,start,end):
    url = 'http://genome.ucsc.edu/cgi-bin/das/hg38/dna?segment={0}:{1},{2}'.format(chr_,start,end)
    response = requests.get(url)
    status_code = response.status_code
    assert status_code == 200
    try:
        my_dict = xmltodict.parse(response.content)
    except:
        print(chr_,start,end)
        raise Exception('Sorry,please check above printed stuff')
    exon_seq = my_dict['DASDNA']['SEQUENCE']['DNA']['#text'].replace('\n','').upper()
    return exon_seq
    
    
    
    
    
    
def fasta_to_dict(path):
    dict_fa = {}
    with open(path,'r') as in_handle:
        for title,seq in SimpleFastaParser(in_handle):
            temp_list = []
            EnsID = title.split('|')[0]
            chro = title.split('|')[1]
            start = title.split('|')[2]
            end = title.split('|')[3]
            temp_list=[chro,start,end,seq]
            dict_fa[EnsID] = temp_list
    return dict_fa
        
def query_from_dict_fa(dict_fa,abs_start,abs_end,EnsID,strand):
    if strand == '+':        
        start = int(dict_fa[EnsID][1])
        end = int(dict_fa[EnsID][2])
        seq = dict_fa[EnsID][3]
        start_index = int(abs_start) - start + 2000
        end_index = int(abs_end) - start + 1 + 2000
        exon_seq = seq[start_index:end_index]
    
    elif strand == '-':
        start = int(dict_fa[EnsID][1])
        end = int(dict_fa[EnsID][2])
        seq_reverse = dict_fa[EnsID][3]
        seq_forward = str(Seq(seq_reverse).reverse_complement())  # Hs_gene.fa restore the reverse strand info
        start_index = int(abs_start) - start + 2000
        end_index = int(abs_end) - start + 1 + 2000 # endpoint in python is non-inclusive
        exon_seq_1 = seq_forward[start_index:end_index]
        s = Seq(exon_seq_1)
        exon_seq = str(s.reverse_complement())
    return exon_seq

def exonCoords_to_dict(path,delimiter):
    coords=[]
    dict_exonCoords={}
    with open(path,'r') as file:
        next(file)
        for line in file:
            items = line.split('\t')
            coords=(items[2],items[3],items[4],items[5])
            if items[0] in dict_exonCoords:
                dict_exonCoords[items[0]][items[1]] = coords
            else:
                dict_exonCoords[items[0]] = {}
                dict_exonCoords[items[0]][items[1]] = coords
    # final structure {'EnsID':{E1:[chr,strand,start,end],E2:[chr,strand,start,end]}}
    return dict_exonCoords

 
def convertExonList(df):
    dictExonList = {}
    for i in range(df.shape[0]):
        EnsGID = df.iat[i,0]
        EnsTID = df.iat[i,1]
        exonList = df.iat[i,3]
        try: dictExonList[EnsGID][EnsTID] = exonList
        except KeyError: dictExonList[EnsGID] = {EnsTID:exonList}
        # {EnsGID:{EnsTID:exonlist,EnsTID:exonlist}}
    return dictExonList
    


    

def spread(list_):
    ret = []
    for i in list_:
        ret.extend(i) if isinstance(i,list) else ret.append(i)
    ret = set(ret)
    return ret

    
def intron(event,EnsGID,junction,dict_exonCoords,dictExonList,N):
    merBucket = []
    if event.startswith('E'):   # only consider this situation, since intron retentions are predominantly associated with a translatable preceding exon, ending up with a early stop codon
        # namely: E2.4-I2.1, E22.1-I22.1
        former = event.split('-')[0]  # E2.4, E22.1
        latter = event.split('-')[1]   # I2.1
        #print(event,EnsGID,former)
        #try: attrs_former = dict_exonCoords[EnsGID][former] # chr, strand, start, end
        #except KeyError: print('preexisting bug, Event {0} in {1} doesn\'t have {2}!!!!!'.format(event,EnsGID,former))
        #else:
        attrs_latter = dict_exonCoords[EnsGID][latter]
        strand = attrs_latter[1]        
        allTransDict = dictExonList[EnsGID] 
        for tran,exonlist in allTransDict.items():
            if former in exonlist: 
                try:tranStartIndex = grabEnsemblTranscriptTable(tran)
                except HTTPError: continue   # for instance, BC4389439 it is not a Ensemble-valid transcript ID, just pass this transcript
                else:
                    if type(tranStartIndex) == int:
                        if strand == '+':
                            intronStartIndex = int(attrs_latter[2])   # index by intron itself
                            remainder = (intronStartIndex - tranStartIndex) % 3   # how many nts remaining before the first nt in intron
                            # if 0, means the first nt in intron will be the first nt in codon triplet,
                            # if 1, means the first nt in intron will be the second nt in codon triplet,
                            # if 2, means the first nt in intron will be the third nt in codon triplet.
                            if remainder == 0: 
                                front=0 if junction.find(',') - ((N-1)*3)<0 else junction.find(',') - ((N-1)*3)
                                junctionSlice = junction[front:].replace(',','')
                            elif remainder == 1: 
                                front=0 if junction.find(',') - ((N-1)*3+1)<0 else junction.find(',') - ((N-1)*3+1)
                                junctionSlice = junction[front:].replace(',','')
                            elif remainder == 2: 
                                front=0 if junction.find(',') - ((N-1)*3+2)<0 else junction.find(',') - ((N-1)*3+2)
                                junctionSlice = junction[front:].replace(',','')
                        elif strand == '-':
                            intronStartIndex = int(attrs_latter[3]) 
                            remainder = (tranStartIndex - intronStartIndex) % 3 
                            if remainder == 0: 
                                front=0 if junction.find(',') - ((N-1)*3)<0 else junction.find(',') - ((N-1)*3)
                                junctionSlice = junction[front:].replace(',','')
                            elif remainder == 1: 
                                front=0 if junction.find(',') - ((N-1)*3+1)<0 else junction.find(',') - ((N-1)*3+1)
                                junctionSlice = junction[front:].replace(',','')
                            elif remainder == 2: 
                                front=0 if junction.find(',') - ((N-1)*3+2)<0 else junction.find(',') - ((N-1)*3+2)
                                junctionSlice = junction[front:].replace(',','')
                
                        try:peptide = str(Seq(junctionSlice).translate(to_stop=False))
                        except:print('**** preexising bugs')
                        else:
                            merArray = extractNmer(peptide,N) 
                            merBucket.append(merArray)
        if merBucket == []: merBucket = [[]]  # means no match for the former subexon.
        
    elif event.startswith('I'): merBucket = [[]]
    return merBucket


# # https://rest.ensembl.org/documentation/info/lookup
# def grabEnsemblTranscriptTable(EnsTID):
#     print(EnsTID,'**********************************************************')
#     server = "https://rest.ensembl.org"
#     ext = "/lookup/id/{0}?expand=1".format(EnsTID)     
#     r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})     
#     try: decoded = r.json()
#     except: 
#         sleep(1)   # if not able to grab information, may be rate limited by Ensembl API, 15 times/second, sleep 1s, then retry
#         try:decoded = r.json()
#         except:
#             print('JSON unknown error')
#         else:
#             try: translationStartIndex = decoded['Translation']['start']
#             except KeyError: print('{0} is not translatable'.format(EnsTID)) # might be invalid transcriptID or they don't have tranStartIndex(don't encode protein)
#             else: return translationStartIndex

#     else:
#         try: translationStartIndex = decoded['Translation']['start']
#         except KeyError: print('{0} is not translatable'.format(EnsTID)) # might be invalid transcriptID or they don't have tranStartIndex(don't encode protein)
#         else: return translationStartIndex
# # for except branch, we don't specify return command, so it will return NoneType   


def dictGTFconstruct(dataFolder):
    dictGTF = {}
    gtfEnsembl91 = pd.read_csv(os.path.join(dataFolder,'gtfEnsembl91.txt'),sep='\t')
    for i in range(gtfEnsembl91.shape[0]):
        feature = gtfEnsembl91.iloc[i]['feature']
        start = gtfEnsembl91.iloc[i]['start']
        tranID = gtfEnsembl91.iloc[i]['tranID']
        if feature == 'start_codon':
            dictGTF[tranID] = start
    return dictGTF


    
def grabEnsemblTranscriptTable(EnsTID):   # replace old grabEnsemblTranscriptTable function with a pre-downloaded dataset
    try:
        translationStartIndex = dictGTF[EnsTID]
    except KeyError:
        print('{0} does not have valid EnsTID'.format(EnsTID))
    else: return translationStartIndex



def toFasta(list_,N):
    with open('./resultMHC/queryNetMHC_{0}.fa'.format(N),'w') as file1:
        for index,item in enumerate(list_):
            file1.write('>{0} mer\n'.format(index+1))
            file1.write('{0}\n'.format(item.strip('\'')))

def inspectGTEx(event,dicTissueExp,tissue='all',plot=True):
    flag = 0
    if tissue=='all':
        tissueExp = dicTissueExp[event]
        for tis,exp in tissueExp.items():
            exp = exp.astype('float64')
            exp=exp[np.logical_not(np.isnan(exp))]
            if exp.size == 0: print('{0} data incomplete'.format(tis))
            elif np.any(exp):   # have non-zero element
                if plot==True:
                    fig = plt.figure()
                    plt.bar(np.arange(len(exp)),exp,width=0.2,label=tis)
                    plt.xticks(np.arange(len(exp)),np.arange(len(exp))+1)
                    plt.legend()
                    plt.savefig(os.path.join(outFolder,'resultMHC_{0}'.format(taskName),'figures/{0}_{1}.pdf'.format(event,tis)),bbox_inches='tight')
                    plt.close(fig)
                else: continue
            else: 
                flag += 1
                print('No expression in {}'.format(tis))
            
    else:
        expression = dicTissueExp[event][tissue]
        exp = expression.astype('float64')
        exp=exp[np.logical_not(np.isnan(exp))]
        if exp.size == 0: print('{0} data incomplete'.format(tissue))
        elif np.any(exp):   # have non-zero element
            plt.bar(np.arange(len(exp)),exp,width=0.2,label='tissue')
            plt.legend()
            plt.savefig(os.path.join(outFolder,'resultMHC_{0}'.format(taskName),'figures/{}.pdf'.format(tissue)),bbox_inches='tight')
            plt.show()
            print(expression)
    return flag  


def run(NeoJBaml):
    PID = os.getpid()

    print('retrieve all junction site sequence-Process:{0}\n'.format(PID))
    NeoJBaml.retrieveJunctionSite(dict_exonCoords,dict_fa)
 
    print('inspecting each event see if they could match up with any existing transcripts-Process:{0}\n'.format(PID))
    NeoJBaml.matchWithExonlist(df_exonlist,dict_exonCoords)
    NeoJBaml.rescueEvent_1()
    NeoJBaml.rescueEvent_2()
 
    print('getting most likely ORF and ORFaa for each event that could match up with existing transcript-Process:{0}\n'.format(PID))
    NeoJBaml.getORF()
    NeoJBaml.getORFaa()

    print('checking the phase of each event\'s junction site-Process:{0}\n'.format(PID))
    NeoJBaml.phaseTranslateJunction()

    print('getting Nmer, only for first round-Process:{0}\n'.format(PID))
    NeoJBaml.getNmer()

    print('first round ends-Process:{0}\n'.format(PID))

    print('starting to mannal check,second round-Process:{0}\n'.format(PID))
    NeoJBaml.mannual()
    print('finished mannual check-Process:{0}\n'.format(PID))

    if mode == 'Neoantigen':
        print('starting to deploy binding affinity prediction-Process:{0}\n'.format(PID))
        if MHCmode == 'MHCI':
            #NeoJBaml.MHCflurry(HLA,MHCmode)
            NeoJBaml.netMHCresult(HLA,software,MHCmode)
        elif MHCmode == 'MHCII':
            NeoJBaml.netMHCresult(HLA,software,MHCmode) 
        print('finished binding affinity prediction-Process:{0}\n'.format(PID))
        return NeoJBaml.df
    elif mode == 'Peptide':
        return NeoJBaml.df   # only get splicing junction peptides
    elif mode == 'immuno+':
        print('starting to deploy binding affinity prediction-Process:{0}\n'.format(PID))
        #NeoJBaml.MHCflurry(HLA,MHCmode)
        NeoJBaml.netMHCresult(HLA,software,MHCmode)
        print('finished binding affinity prediction-Process:{0}\n'.format(PID))
        print('starting immunogenecity prediction\n')
        NeoJBaml.immunogenecity(HLA)
        print('finishing immunogenecity prediction\n')
        return NeoJBaml.df




def main(intFile,taskName,outFolder,dataFolder,k,HLA,software,MHCmode,mode,Core,checkGTEx,cutoff_PSI,cutoff_sampleRatio,cutoff_tissueRatio):

    startTime = process_time()
    global df
    global df_exonlist
    global dict_exonCoords
    global dict_fa
    global dictExonList
    global dictGTF
    global hla
    global dic_inventory
    #global predictor
    if not os.path.exists(os.path.join(outFolder,'resultMHC_{0}'.format(taskName))): os.makedirs(os.path.join(outFolder,'resultMHC_{0}'.format(taskName)))
    if not os.path.exists(os.path.join(outFolder,'resultMHC_{0}'.format(taskName),'temp')): os.makedirs(os.path.join(outFolder,'resultMHC_{0}'.format(taskName),'temp'))
    if not os.path.exists(os.path.join(outFolder,'resultMHC_{0}'.format(taskName),'figures')): os.makedirs(os.path.join(outFolder,'resultMHC_{0}'.format(taskName),'figures'))
    print('loading input files\n')
    df = pd.read_csv(intFile,sep='\t') 
    print('loading all existing transcripts files\n')
    df_exonlist = pd.read_csv(os.path.join(dataFolder,'mRNA-ExonIDs.txt'),sep='\t',header=None,names=['EnsGID','EnsTID','EnsPID','Exons'])
    print('loading all subexon coordinates files\n')
    dict_exonCoords = exonCoords_to_dict(os.path.join(dataFolder,'Hs_Ensembl_exon.txt'),'\t')
    print('loading exon sequence fasta files, 2GB\n')
    dict_fa = fasta_to_dict(os.path.join(dataFolder,'Hs_gene-seq-2000_flank.fa'))
    print('converting subexon coordinates to a dictionary\n')
    dictExonList = convertExonList(df_exonlist)
    print('constructing dictGTF file')
    dictGTF = dictGTFconstruct(dataFolder)
    print('some necessary files for immunogenecity prediction\n')
    # hla = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/immunogenecity/transformer/hla2paratopeTable_aligned.txt',sep='\t',header=None,names=['hla','paratope'])
    # inventory = hla['hla']
    # dic_inventory = dict_inventory(inventory)


  
    if checkGTEx == 'True':
        if os.path.exists(os.path.join(outFolder,'resultMHC_{0}'.format(taskName),'{0}_neojunction_singleSample_check.txt'.format(taskName))):
            dfNeoJunction = pd.read_csv(os.path.join(outFolder,'resultMHC_{0}'.format(taskName),'{0}_neojunction_singleSample_check.txt'.format(taskName)),sep='\t')
            print('NeoJunctions are ready\n')
        else:
            global dicTissueExp
            print('loading GTEx dataset, it will take 20 mins, please be patient')
            start = process_time()
            with bz2.BZ2File(os.path.join(dataFolder,'dicTissueExp.pbz2'),'rb') as f1:
                dicTissueExp = cpickle.load(f1)  
                end = process_time()    
            print('consume {0}'.format(end-start))
            metaBaml = Meta(df) #Instantiate Meta object
            print('generate NeoJunctions\n')
            dfNeoJunction = check_GTEx_PSI_pbz2(metaBaml.df,cutoff_PSI,cutoff_sampleRatio,cutoff_tissueRatio,taskName,outFolder,dataFolder)
            if dfNeoJunction.shape[0] == 0:
                raise Exception('After checking GTEx, no event remains')
            dfNeoJunction.to_csv(os.path.join(outFolder,'resultMHC_{0}'.format(taskName),'{0}_neojunction_singleSample_check.txt'.format(taskName)),sep='\t',index=None)
            print('NeoJunctions are ready\n')
        
    if checkGTEx == 'False':
        metaBaml = Meta(df) #Instantiate Meta object
        print('generate NeoJunctions\n')
        dfNeoJunction = neoJunction_testMethod_noGTEx(metaBaml.df)
        print('NeoJunctions are ready\n')

    if mode == 'GTEx_check':
        sys.exit('Congrats, Already successfully finished, check the output files for GTEx check and IW scores.')
       
    print('start analysis and spawn subprocesses\n')




    df_split = np.array_split(dfNeoJunction, Core, axis=0)    # cut by rows, equally splited 
    obj_split = [NeoJ(df,k) for df in df_split]
    print('start pooling\n')

    pool = Pool(Core)
    df_out_list = pool.map(run, obj_split)
    pool.close()
    pool.join()
    Crystal = pd.concat(df_out_list)   # default is stacking by rows
    Crystal = Crystal.drop(['exam_seq','exam_whole_transcripts','exam_ORF_tran','exam_ORF_aa','phase'],axis=1)
    if mode == 'Peptide':
        Crystal.to_csv(os.path.join(outFolder,'resultMHC_{1}/JunctionPeptides_{0}_{1}.txt'.format(k,taskName)),sep='\t',header=True,index = False)
    elif mode == 'Neoantigen':
        Crystal.to_csv(os.path.join(outFolder,'resultMHC_{1}/Neoantigen_{0}_{1}.txt'.format(k,taskName)),sep='\t',header=True,index = False)
    elif mode == 'immuno+':
        Crystal.to_csv(os.path.join(outFolder,'resultMHC_{1}/immuno+_{0}_{1}.txt'.format(k,taskName)),sep='\t',header=True,index = False)

    
    endTime = process_time()
    print('Time Usage: {} seconds'.format(endTime-startTime))




def usage():


    print('Usage:')
    print('python3 mhcPresent.py -i /data/salomonis2/LabFiles/Frank-Li/breast_cancer_run/all_events/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt -t breast_all -o /data/salomonis2/LabFiles/Frank-Li/breast_cancer_run/all_events -d /data/salomonis2/LabFiles/Frank-Li/python3/data -k 8 -H HLA-A01:01,HLA-A03:01,HLA-B07:02,HLA-B27:05,HLA-B58:01 -s /data/salomonis2/LabFiles/Frank-Li/python3/netMHCpan-4.1/netMHCpan -M MHCI -m Peptide -C 8 -c False --cutoffPSI 0.1 --cutoffSample 0.1 --cutoffTissue 0.1')
    print('Options:')
    print('-i : path of input file')
    print('-t : the name of your task')
    print('-o : output folder')
    print('-d : data folder')
    print('-k : Kmer')
    print('-H : queried HLA alleles')
    print('-s : full path for where netMHCpan sit')
    print('-M : MHCmode, either MHCI or MHCII')
    print('-m : Get Splicing peptides or Neoantigens')
    print('-C : how many processes you wanna spawn')
    print('-c: check GTEx data or not')
    print('--cutoffPSI : above this PSI value will be labelled as expressed')
    print('--cutoffSample : above this ratio this tissue will be labelled as expressed')
    print('--cutoffTissue : above this ratio this event will be labelled as expressed in normal tissue in general')
    print('-h --help : check help information ')
    print('Author: Guangyuan(Frank) Li <li2g2@mail.uc.edu>, PhD Student, University of Cincinnati, 2020') 
    



if __name__ == "__main__":

#    log_err = open('queryGTEx.stderr.log','a')
#    log_out = open('queryGTEx.stdout.log','a')
#    sys.stderr = log_err
#    sys.stdout = log_out
    try:
        options, remainder = getopt.getopt(sys.argv[1:],'hi:t:o:d:k:H:s:M:m:C:c:',['help','cutoffPSI=','cutoffSample=','cutoffTissue='])
    except getopt.GetoptError as err:
        print('ERROR:', err)
        usage()
        sys.exit(1)
    for opt, arg in options:
        if opt in ('-i'):
            intFile = str(arg)
            print('Input file is:', arg)
        elif opt in ('-t'):
            taskName = str(arg)
            print('give your task a name:',arg)
        elif opt in ('-o'):
            outFolder = str(arg)
            print('output folder:',arg)
        elif opt in ('-d'):
            dataFolder = str(arg)
            print('data folder:',arg)
        elif opt in ('-k'):
            k= int(arg)
            print('kmer:', arg)
        elif opt in ('-H'):
            HLA = str(arg)
            print('Queried HLA allele:',arg)
        elif opt in ('-s'):
            software = str(arg)
            print('Full path of netMHCpan standalone version:',arg)
        elif opt in ('-M'):
            MHCmode = str(arg)
            print('MHCI or MHCII:',arg)
        elif opt in ('-m'):
            mode = str(arg)
            print('Get Splicing peptides or Neoantigens:',arg)
        elif opt in ('-C'):
            Core = int(arg)
            print('How many processes to use:',arg)
        elif opt in ('-c'):
            checkGTEx = arg
            print('check GTEx or not:',arg)
        elif opt in ('--cutoffPSI'):
            cutoff_PSI = float(arg)
            print('cutoff for PSI value:',arg)
        elif opt in ('--cutoffSample'):
            cutoff_sample = float(arg)
            print('cutoff for sample ratio:',arg)    
        elif opt in ('--cutoffTissue'):
            cutoff_tissue = float(arg)
            print('cutoff for tissue ratio:',arg) 
        elif opt in ('--help','-h'):
            usage() 
            sys.exit() 


    main(intFile,taskName,outFolder,dataFolder,k,HLA,software,MHCmode,mode,Core,checkGTEx,cutoff_PSI,cutoff_sample,cutoff_tissue)

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    




