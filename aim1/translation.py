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

# other modules
from binding import run_netMHCpan
from deepimmuno import run_deepimmuno

'''
Now for a junction, you need to obtain the translated neo-epitopes, simply put, you just need two things
1. junction DNA sequence
2. How to translate them
'''

class EnhancedPeptides():
    def __init__(self,peptides,hlas,kind):
        self.mers = []
        self.info = []
        for k,v in peptides.items():
            self.mers.append(k)
            phlas = {}  # {pep1:{hla1:{}}}
            for pep in v:
                if kind == 0:  # each pep is just a string
                    pairs = {}
                    for hla in hlas:
                        pairs[hla] = {}
                    phlas[pep] = pairs
                elif kind == 1:   # echo pep is (pep,hla)
                    try: 
                        phlas[pep[0]][pep[1]] = {}
                    except KeyError:
                        phlas[pep[0]] = {}
                        phlas[pep[0]][pep[1]] = {}
            self.info.append(phlas)

    def __str__(self):
        return str(self.info)

    def register_attr(self,df,attr_name):
        '''
        df should follow:
          peptide     mer     hla       score    identity
         AAAAAAAAA     9   HLA-A*01:01   0.3       SB        
        '''
        for mer,sub_df in df.groupby(by='mer'):
            index = self.mers.index(int(mer))
            for peptide,sub_sub_df in sub_df.groupby(by='peptide'):
                for row in sub_sub_df.itertuples(index=False):
                    self.info[index][peptide][row.hla][attr_name] = (row.score,row.identity)

    def filter_based_on_criterion(self,criteria):
        # criterion: [(net),], ['netMHCpan_el',1,==,SB]
        peptides = {k:[] for k in self.mers}
        for i,k in enumerate(self.mers):
            for pep,hla_complex in self.info[i].items():
                for hla,attrs in hla_complex.items():
                    boolean_list = []
                    for criterion in criteria:
                        boolean = eval('attrs[\'{}\'][{}] {} \'{}\''.format(criterion[0],criterion[1],criterion[2],criterion[3]))
                        boolean_list.append(boolean)
                    boolean_final = all(boolean_list)
                    if boolean_final:
                        peptides[k].append((pep,hla))
        return peptides




        


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

    def in_silico_translation(self,ks=[9,10]):
        peptides = {k:[] for k in ks}
        if '$' not in self.junction and '*' not in self.junction and '#' not in self.junction:
            first,second = self.junction.split(',')
            for phase in [0,1,2]:  # tranlation starts with index "phase"
                de_facto_first = first[phase:]
                pep_dict = get_peptides(de_facto_first,second,ks)
                for k in ks:
                    peptides[k].extend(pep_dict[k])
                    peptides[k] = list(set(peptides[k]))
        self.peptides = peptides
        return peptides

    def binding_prediction(self,hlas):
        ep = EnhancedPeptides(self.peptides,hlas,0)
        for k,v in self.peptides.items():
            df = run_netMHCpan(software_path,v,hla_formatting(hlas,'netMHCpan_output','netMHCpan_input'),k)
            ep.register_attr(df,'netMHCpan_el')
        self.enhanced_peptides = ep
        return ep

    def immunogenecity_prediction(self):
        reduced = self.enhanced_peptides.filter_based_on_criterion([('netMHCpan_el',1,'==','SB'),])
        ep = EnhancedPeptides(reduced,hlas,1)
        for k,v in reduced.items():
            data = np.array(v)
            df_input = pd.DataFrame(data=data)
            df_input[1] = hla_formatting(df_input[1].tolist(),'netMHCpan_output','deepimmuno')
            df_output = run_deepimmuno(df_input)
            df = pd.DataFrame({'peptide':df_output['peptide'].values,'mer':[k]*df_output.shape[0],
                               'hla':hla_formatting(df_output['HLA'].values.tolist(),'deepimmuno','netMHCpan_output'),
                               'score':df_output['immunogenicity'].values,
                               'identity':[True if item > 0.5 else False for item in df_output['immunogenicity'].values]})
            ep.register_attr(df,attr_name='deepimmuno_immunogenicity')
            self.enhanced_peptides.register_attr(df,attr_name='deepimmuno_immunogenicity')
        return ep
        
            

# processing functions
def hla_formatting(pre,pre_type,post_type):
    if pre_type == 'netMHCpan_output' and post_type == 'netMHCpan_input':  # HLA-A*01:01 to HLA-A01:01
        post = [hla.replace('*','') for hla in pre]
    elif pre_type == 'netMHCpan_output' and post_type == 'deepimmuno':  # HLA-A*01:01 to HLA-A*0101
        post = [hla.replace(':','') for hla in pre]
    elif pre_type == 'deepimmuno' and post_type == 'netMHCpan_output': # HLA-A*0101 to HLA-A*01:01
        post = [hla[:8] + ':' + hla[-2:] for hla in pre]
    return post





def get_peptides(de_facto_first,second,ks):
    peptides = {k:[] for k in ks}
    extra = len(de_facto_first) % 3  # how many base left in first assuming no stop condon in front of it.
    num = len(de_facto_first) // 3   # how many set of codons in the first.
    aa_first = str(Seq(de_facto_first).translate(to_stop=True))
    if len(aa_first) == num:  # successfully read through
        if extra == 0:
            continue_second = second
        elif extra == 1:
            continue_second = de_facto_first[-1] + second
        elif extra == 2:
            continue_second = de_facto_first[-2:] + second
        aa_second = str(Seq(continue_second).translate(to_stop=True))
        if len(aa_second) > 0:  # at least, not ''
            for k in ks:
                second_most = min(k,len(aa_second)) # the max allowed number of base for aa_second
                first_most = len(aa_first)  # the max allowed number of base for aa_first
                for n_from_second in range(second_most,0,-1):
                    n_from_first = k - n_from_second
                    if n_from_first == 0 and extra == 0:
                        '''
                        cttca cct cac ttt acc ttc tcc tcc agc aca gga act agg aac tac gga gag aga agc caa 
                           S   P   H   F   T   F   S   S   S   T   G   T   R   N   Y   G   E   R   S   Q

                        the common is between "ttt" and "acc"
                        so, when extra == 0, means whole second will be after the junction, in this case, we have to require 
                        the peptide at least include a amino acid from first part, otherwise, TFSSSTGTR won't be a splice-derived
                        peptide.   

                        this example can be reproduced using: nj = NeoJunction('ENSG00000223572:E15.1-E15.2')                     
                        '''
                        continue
                    if n_from_first <= first_most:
                        if n_from_first > 0:
                            pep = aa_first[-n_from_first:] + aa_second[:n_from_second]
                        elif n_from_first == 0:
                            pep = aa_second[:n_from_second]
                        peptides[k].append(pep)
    return peptides
                        


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
    software_path = '../external/netMHCpan-4.1/netMHCpan'
    hlas = ['HLA-A*01:01','HLA-A*02:01','HLA-A*24:02','HLA-A*68:01','HLA-B*08:01','HLA-B*08:02']
    nj = NeoJunction('ENSG00000223572:E15.1-E15.2')
    nj.detect_type()
    nj.retrieve_junction_seq()
    nj.in_silico_translation()
    nj.binding_prediction(hlas=hlas)
    nj.immunogenecity_prediction()




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
    





