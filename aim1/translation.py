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
import re
import requests
import xmltodict

# for biopython, pip install biopython
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq

'''
Now for a junction, you need to obtain the translated neo-epitopes, simply put, you just need two things
1. junction DNA sequence
2. How to translate them
'''

class NeoJunction():
    def __init__(self,uid):
        self.uid = uid
    
    def detect_type(self):
        '''
        Ordinary: ENSG00000107902:E10.1-E12.1
        Alt3: ENSG00000110057:E5.1-E6.2_67996641
        Alt5: ENSG00000100321:E7.1_39364266-E8.1
        Intron Retention: ENSG00000115524:I4.1-E5.1
        Novel Exon: ENSG00000008441:I40.1_13076665-E41.1
        Trans-splicing: ENSG00000196565:E14.2-ENSG00000213934:E3.1
        UTR Event: ENSG00000164068:U0.1_49689185-E2.1
        '''
        valid_pattern = re.compile(r'^ENSG\d+:.+?-.+')
        if re.search(valid_pattern,self.uid):   # at least valid one
            if len(re.findall('ENSG',self.uid)) == 2:
                event_type = 'trans_splicing'
            elif 'U' in self.uid:
                event_type = 'utr_event'
            elif '_' in self.uid:
                subexon12 = self.uid.split(':')[1]
                subexon1, subexon2 = subexon12.split('-')
                if 'I' in subexon12:
                    event_type = 'novel_exon'
                elif '_' in subexon1 and '_' in subexon2:
                    event_type = 'alt5_alt3'
                elif '_' in subexon1 and '_' not in subexon2:
                    event_type = 'alt5'
                elif '_' in subexon2 and '_' not in subexon1:
                    event_type = 'alt3'
                else:
                    event_type = 'invalid'
            elif 'I' in self.uid:
                event_type = 'intron_retention'
            elif re.search(r'^ENSG\d+:E\d+\.\d+-E\d+\.\d+$',self.uid):
                event_type = 'ordinary'
            else:
                event_type = 'invalid'
        else:
            event_type = 'invalid'
        self.event_type = event_type
        return event_type

    def retrieve_junction_seq(self):
        if self.event_type != 'invalid':
            ensid = self.uid.split(':')[0]
            subexon1,subexon2 = self.uid.split(':')[1].split('-')
            seq1 = subexon_tran(subexon1,ensid,'site1')
            seq2 = subexon_tran(subexon2,ensid,'site2')
            junction = ','.join([seq1,seq2])
            self.junction = junction
        else:
            self.junction = '$' * 10   # indicating invalid uid
        return junction


# processing functions
def fasta_to_dict(path):
    '''
    Let's talk about the fasta file
    >ENSG|chro|start|end
    seq

    the start and end always correspond to xth base in forward strand,
    however, if the ENSG falls into backward strand, the seq it stored is actually the 
    backward strand from its own 5' - 3'.
    '''
    dict_fa = {}  # {ENSID, [chro,start,end,seq]}
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

def query_from_dict_fa(abs_start,abs_end,EnsID,strand):
    '''
    abs_start and abs_end always means the xth base in forward strand

    the returned exon_seq, however, means the 5'-3' seq depending on the strand information.
    '''
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


'''
Let's just pre-generate the exonCoords with added col for other uses
'''
def is_suffer_from_overhang(path):
    df = pd.read_csv(path,sep='\t')
    store_chunk = []
    for gene,sub_df in df.groupby(by='gene'):
        exonlist = sub_df['exon-id'].tolist()
        dic = {} # {E14 : 1,2,4,5}, the values are stored in str
        for subexon in exonlist:  # this traversal is to build the dict
            exon_num = subexon.split('.')[0]
            subexon_num = subexon.split('.')[1]
            if exon_num in dic:
                dic[exon_num].append(subexon_num)
            else:
                dic[exon_num] = []
                dic[exon_num].append(subexon_num)  
        added_col = []
        for subexon in exonlist:  # this traversal is to figure out whether they are middle subexon or not
            exon_num = subexon.split('.')[0]
            subexon_num = subexon.split('.')[1]
            if str(int(subexon_num) + 1) in dic[exon_num]:
                is_suffer = True
            else:
                is_suffer = False
            added_col.append(is_suffer)
        sub_df['is_suffer'] = added_col
        store_chunk.append(sub_df)
    df = pd.concat(store_chunk,axis=0).sort_index()
    return df
            




def exonCoords_to_dict(path):
    '''
    1. the start and end always forward strand
    2. to clarify the overhang issue, the issue is every middle subexon, its end coord need to be backtracted by 1.
    However, this is different operationally in + and - strand. positive strand is to substract the end by 1 (for middle subexon).
    the negative strand is to add the start by 1 (for middle subexon.)
    '''
    coords=[]
    dict_exonCoords={} # {'EnsID':{E1.1:[chr,strand,start,end,suffer]}} 
    with open(path,'r') as file:
        next(file)
        for line in file:
            items = line.split('\t')
            coords=(items[2],items[3],items[4],items[5],items[10].rstrip('\n'))
            if items[0] in dict_exonCoords:
                dict_exonCoords[items[0]][items[1]] = coords
            else:
                dict_exonCoords[items[0]] = {}
                dict_exonCoords[items[0]][items[1]] = coords
    return dict_exonCoords

def utrAttrs(EnsID):  # try to get U0.1's attribute, but dict_exonCoords doesn't have, so we just wanna get the first entry for its EnsGID
    exonDict = dict_exonCoords[EnsID] 
    attrs = next(iter(exonDict.values()))
    chr_,strand = attrs[0],attrs[1]
    return chr_,strand

def utrJunction(site,EnsGID,strand,chr_,flag,seq_len=100):  # U0.1_438493849, here 438493849 means the site (suffix)
    if flag == 'site1' and strand == '+':  # U0.1_438493849 - E2.2
        otherSite = int(site) - seq_len + 1   # extract UTR with length = 100
        exon_seq = retrieveSeqFromUCSCapi(chr_,int(otherSite),int(site))
    elif flag == 'site1' and strand == '-':    
        otherSite = int(site) + seq_len - 1 
        exon_seq = retrieveSeqFromUCSCapi(chr_,int(site),int(otherSite))
        exon_seq = str(Seq(exon_seq).reverse_complement())
    elif flag == 'site2' and strand == '+':  # E5.3 - U5.4_48374838
        otherSite = int(site) + seq_len -1
        exon_seq = retrieveSeqFromUCSCapi(chr_,int(site),int(otherSite))
    elif flag == 'site2' and strand == '-':
        otherSite = int(site) - seq_len + 1
        exon_seq = retrieveSeqFromUCSCapi(chr_,int(otherSite),int(site))
        exon_seq = str(Seq(exon_seq).reverse_complement())
    return exon_seq

def retrieveSeqFromUCSCapi(chr_,start,end):
    url = 'http://genome.ucsc.edu/cgi-bin/das/hg38/dna?segment={0}:{1},{2}'.format(chr_,start,end)
    response = requests.get(url)
    status_code = response.status_code
    assert status_code == 200
    try:
        my_dict = xmltodict.parse(response.content)
    except:
        exon_seq = '#' * 10  # indicating the UCSC doesn't work
    exon_seq = my_dict['DASDNA']['SEQUENCE']['DNA']['#text'].replace('\n','').upper()
    return exon_seq

def subexon_tran(subexon,EnsID,flag):  # flag either site1 or site2
    '''
    1. subexon can take multiple forms depending on the event type
    E1.2 or I3.4
    E6.2_67996641 or I40.1_13076665, also depending on whether they are subexon1 or subexon2
    ENSG00000213934:E3.1 or ENSG00000213934:E2.1_473843893894
    U0.1_49689185

    2. everything with trailing suffix will depend on the subexon1 or subexon2, but sometimes, it is fixed (trans-splicing can only be in subexon2)
    3. to be clear, the exon_seq returned is always 5'-3' sequence, not forward anymore.

    '''
    try:   # E1.2 or I3.4
        attrs = dict_exonCoords[EnsID][subexon]  # [chr,strand,start,end,suffer]
        if attrs[1] == '+':  
            if attrs[4] == 'True':  # remedy by substract the end by 1
                exon_seq = query_from_dict_fa(attrs[2],str(int(attrs[3])-1),EnsID,attrs[1]) 
            else:
                exon_seq = query_from_dict_fa(attrs[2],attrs[3],EnsID,attrs[1]) 
        else:   
            if attrs[4] == 'True': # remedy by adding the start by 1
                exon_seq = query_from_dict_fa(str(int(attrs[2])+1),attrs[3],EnsID,attrs[1]) 
            else:
                exon_seq = query_from_dict_fa(attrs[2],attrs[3],EnsID,attrs[1]) 
    except KeyError:
        if ':' in subexon: # ENSG00000213934:E3.1
            fusionGeneEnsID = subexon.split(':')[0] 
            fusionGeneExon = subexon.split(':')[1]        
            if  '_' in fusionGeneExon:   # ENSG:E2.1_473843893894
                suffix = fusionGeneExon.split('_')[1]
                subexon = fusionGeneExon.split('_')[0]
                attrs = dict_exonCoords[fusionGeneEnsID][subexon]
                if attrs[1] == '+':  
                    if attrs[4] == 'True': # remedy by substracting the end by 1
                        exon_seq = query_from_dict_fa(suffix,str(int(attrs[3])-1),fusionGeneEnsID,attrs[1]) 
                    else:
                        exon_seq = query_from_dict_fa(suffix,attrs[3],fusionGeneEnsID,attrs[1]) 
                else:  
                    if attrs[4] == 'True':  # remedy by adding the start by 1
                        exon_seq = query_from_dict_fa(str(int(attrs[2])+1),suffix,fusionGeneEnsID,attrs[1])
                    else:
                        exon_seq = query_from_dict_fa(attrs[2],suffix,fusionGeneEnsID,attrs[1])
            else:  # ENSG:E2.1
                try:
                    attrs = dict_exonCoords[fusionGeneEnsID][fusionGeneExon]
                except KeyError:
                    exon_seq = '*' * 10  # indicator for error on MultiPath-PSI itself
                else:
                    if attrs[1] == '+':  
                        if attrs[4] == 'True':  # remedy by substract the end by 1
                            exon_seq = query_from_dict_fa(attrs[2],str(int(attrs[3])-1),fusionGeneEnsID,attrs[1]) 
                        else:
                            exon_seq = query_from_dict_fa(attrs[2],attrs[3],fusionGeneEnsID,attrs[1]) 
                    else:   
                        if attrs[4] == 'True': # remedy by adding the start by 1
                            exon_seq = query_from_dict_fa(str(int(attrs[2])+1),attrs[3],fusionGeneEnsID,attrs[1]) 
                        else:
                            exon_seq = query_from_dict_fa(attrs[2],attrs[3],fusionGeneEnsID,attrs[1]) 

        else:  # could be trailing or utr, or non-existing ordinary subexon
            try:
                suffix = subexon.split('_')[1]
            except IndexError: # the logic is there's a subexon E45.3, it is no trailing, but just not in the exonCoords.
                exon_seq = '*' * 10  # indicator for error on MultiPath-PSI itself
            else:
                subexon = subexon.split('_')[0]
                try:
                    attrs = dict_exonCoords[EnsID][subexon]
                except KeyError:  # must be UTR
                    chrUTR,strandUTR = utrAttrs(EnsID) # this is get from a random subexon under that EnsID
                    exon_seq = utrJunction(suffix,EnsID,strandUTR,chrUTR,flag)  
                else:   # must be trailing
                    if flag == 'site2':
                        if attrs[1] == '+':  
                            if attrs[4] == 'True': # remedy by substracting the end by 1
                                exon_seq = query_from_dict_fa(suffix,str(int(attrs[3])-1),EnsID,attrs[1]) 
                            else:
                                exon_seq = query_from_dict_fa(suffix,attrs[3],EnsID,attrs[1]) 
                        else:  
                            if attrs[4] == 'True':  # remedy by adding the start by 1
                                exon_seq = query_from_dict_fa(str(int(attrs[2])+1),suffix,EnsID,attrs[1])
                            else:
                                exon_seq = query_from_dict_fa(attrs[2],suffix,EnsID,attrs[1])
                    elif flag == 'site1':  # not affected by overhang since it is site1
                        if attrs[1] == '+': 
                            exon_seq = query_from_dict_fa(attrs[2],suffix,EnsID,attrs[1])
                        else:
                            exon_seq = query_from_dict_fa(suffix,attrs[3],EnsID,attrs[1])
    return exon_seq




if __name__ == '__main__':
    dict_exonCoords = exonCoords_to_dict('../data/Hs_Ensembl_exon_add_col.txt')
    dict_fa = fasta_to_dict('../data/Hs_gene-seq-2000_flank.fa')
    exon_seq = subexon_tran('E2.1_105625354','ENSG00000211892','site1')
    print(exon_seq)


    sys.exit('stop')
    df = is_suffer_from_overhang('../data/Hs_Ensembl_exon.txt')
    df.to_csv('../data/Hs_Ensembl_exon_add_col.txt',sep='\t',index=None)





    # testing case 
    # some quick tips
    '''
    when xth, the direct subtraction means from one end, how many to count to arrive the another end
    0 based index, always means how many base before it.
    '''
    
    # ordinary + 
    subexon_tran('E15.1','ENSG00000223572','site1')
    subexon_tran('E15.2','ENSG00000223572','site1')

    # ordinary -
    subexon_tran('E1.1','ENSG00000149806','site1')
    subexon_tran('E1.2','ENSG00000149806','site1')
    subexon_tran('E1.3','ENSG00000149806','site1')

    # trans with trailing
    subexon_tran('ENSG00000253755:E3.1_105668978','ENSG00000211896','site2')

    # trans without trailing
    subexon_tran('ENSG00000211896:E4.1','ENSG00000211897','site2')

    # utr
    subexon_tran('U0.1_105857815','ENSG00000211895','site1')

    # alt3
    subexon_tran('E11.1_91157640','ENSG00000133943','site2')

    # alt5
    subexon_tran('E2.1_105625354','ENSG00000211892','site1')
    





