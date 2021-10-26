#!/Users/ligk2e/opt/anaconda3/envs/python3/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 19:43:43 2020

@author: ligk2e
"""
import sys
#sys.path.append('/Users/ligk2e/opt/anaconda3/envs/python3/lib/python3.7/site-packages')
import os
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import pandas as pd
from decimal import Decimal as D
import pickle
import bz2
import _pickle as cpickle
import regex
import re
from time import process_time
import collections
import matplotlib.pyplot as plt
import math
import urllib.parse
import urllib.request
import requests
import xmltodict
import subprocess
import argparse
import getopt
import ast
import bisect
import numpy as np






def GetIncreasedPart(df):
    df['sign'] = df['dPSI'].apply(lambda x: True if x>0 else False) # x = lambda a,b:a+b; x(5,6)   
    df_ori = df[df['sign']==True]
    df_ori = df_ori.drop(columns=['sign'])  # how to drop a column
    df_ori = df_ori.set_index(pd.Index(np.arange(df_ori.shape[0])))
    return df_ori




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
        Exons1 = '|' + Exons
        Exons2 = Exons + '|'
        
        if re.search(rf'{re.escape(Exons1)}',item) or re.search(rf'{re.escape(Exons2)}',item) or re.search(rf'{re.escape(Exons)}$',item):
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
                    # expression of minus strand, need to draw an illustrator to visulize that.
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
    #print(strand)
    if strand == '+':
        
        start = int(dict_fa[EnsID][1])
        end = int(dict_fa[EnsID][2])
        seq = dict_fa[EnsID][3]
        start_index = int(abs_start) - start + 2000
        #print(type(start_index))
        end_index = int(abs_end) - start + 1 + 2000
        exon_seq = seq[start_index:end_index]
    
    elif strand == '-':
        start = int(dict_fa[EnsID][1])
        end = int(dict_fa[EnsID][2])
        seq_reverse = dict_fa[EnsID][3]
        seq_forward = str(Seq(seq_reverse,generic_dna).reverse_complement())  # Hs_gene.fa restore the reverse strand info
        start_index = int(abs_start) - start + 2000
        #print(type(start_index))
        end_index = int(abs_end) - start + 1 + 2000 # python range/slice doesn't include end point
        exon_seq_1 = seq_forward[start_index:end_index]
        s = Seq(exon_seq_1,generic_dna)
        exon_seq = str(s.reverse_complement())
    return exon_seq

##############################################################################################    
# part2: pick_peptide.py, find the most likely ORF

###############################################################################################    
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
                 


def exonCoords_to_dict(path,delimiter):
    coords=[]
    dict_exonCoords={}
    with open(path,'r') as file:
        next(file)
        for line in file:
#            dict_temp={}
            items = line.split('\t')
            coords=(items[2],items[3],items[4],items[5])
            if items[0] in dict_exonCoords:
                dict_exonCoords[items[0]][items[1]] = coords
            else:
                dict_exonCoords[items[0]] = {}
                dict_exonCoords[items[0]][items[1]] = coords
    # final structure {'EnsID':{E1:[chr,strand,start,end],E2:[chr,strand,start,end]}}
    return dict_exonCoords

            

#################################################################################################
# part3: following.py   find the junction sites' sequence.   
#################################################################################################

def retrieveJunctionSite(df,dict_exonCoords,dict_fa):
    exam_seq,back_seq = [],[]
    for i in range(df.shape[0]):
        temp = uid(df,i)
        EnsID = list(temp.keys())[0].split(':')[1]
        exam_site = list(temp.values())[0][0]
        back_site = list(temp.values())[0][1]
        exam_site_1 = subexon_tran(exam_site.split('-')[0],EnsID,dict_exonCoords,dict_fa,'site1')
        exam_site_2 = subexon_tran(exam_site.split('-')[1],EnsID,dict_exonCoords,dict_fa,'site2')
        exam_seq_join = ','.join([exam_site_1,exam_site_2])
        exam_seq.append(exam_seq_join)
        back_site_1 = subexon_tran(back_site.split('-')[0],EnsID,dict_exonCoords,dict_fa,'site1')
        back_site_2 = subexon_tran(back_site.split('-')[1],EnsID,dict_exonCoords,dict_fa,'site2')
        back_seq_join = ','.join([back_site_1,back_site_2])
        back_seq.append(back_seq_join)
        
    df['exam_seq'] = exam_seq
    #df['back_seq'] = back_seq
    return df
        
        
def subexon_tran(subexon,EnsID,dict_exonCoords,dict_fa,flag): # flag means if it is site_1 or site_2
    try:
        attrs = dict_exonCoords[EnsID][subexon]
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
            except IndexError:
                exon_seq = '***********************'   
                print('{0} does not include in {1} exonlists'.format(subexon,EnsID))
            else:
                subexon = subexon.split('_')[0]
                try:
                    attrs = dict_exonCoords[EnsID][subexon]
                except KeyError:
                    print('{0} observes an UTR event {1}'.format(EnsID,subexon))
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
        exon_seq = str(Seq(exon_seq,generic_dna).reverse_complement())
    elif flag == 'site2' and strand == '+':  # E5.3 - U5.4_48374838
        otherSite = int(site) + 100 -1
        exon_seq = retrieveSeqFromUCSCapi(chr_,int(site),int(otherSite))
    elif flag == 'site2' and strand == '-':
        otherSite = int(site) - 100 + 1
        #print(EnsGID,chr_,site,otherSite)
        exon_seq = retrieveSeqFromUCSCapi(chr_,int(otherSite),int(site))
        exon_seq = str(Seq(exon_seq,generic_dna).reverse_complement())
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


##################################################################################################
# part4: narrow down to extracellular instances and check if good representative and do seeding alignment   
###################################################################################################
    

 


                        
def check_if_good_representative(df):
    outer_condition_array = []
    outer_condition_detail = []
    
    
    for i in range(df.shape[0]):
        condition_array = []
        first = df.iloc[i]['exam_first_whole_transcripts']
        second = df.iloc[i]['second_round']
        third =df.iloc[i]['third_round']
        if list_check(first): all_whole_tran = first
        elif not list_check(first) and list_check(second): all_whole_tran = second
        elif not list_check(first) and not list_check(second) and list_check(third): all_whole_tran = third
        else:
            all_whole_tran = 'grenade'
            condition_array=[False]   # intron, unrecoverable or third round is still all empty, will fall into df_filtered
      

        if not all_whole_tran == 'grenade':



            all_ORF_tran = df['ORF'].tolist()[i]
            junction = df['exam_seq'].tolist()[i].replace(',','')

            for idx in range(len(all_whole_tran)):
                whole = all_whole_tran[idx]
                ORF = all_ORF_tran[idx]
                if whole:
                    whole = whole.replace(',','')
                    start_ORF = whole.find(ORF)
                    end_ORF = start_ORF + len(ORF)
                    start_junction = whole
        # we have to apply fuzzy matching, because junction consists of former part and latter part(i.e. E6.3|E8.1)
        # In my way to handle overlapping 1nt, the former one will always stay constant, but latter one in the whole
        # transcript, it might get trimmed but here we don't trim it, so there might be 1 overhang in jucntion seq.
        
                    pattern = regex.compile('(%s){d<=1}' % junction) 
                    try:
                        start_junction = pattern.search(whole).span()[0]
                    except:
                        print(df.iloc[i]['UID'],idx)
                        raise Exception('bug')
                    end_junction = pattern.search(whole).span()[1] - 1

                    if start_junction <= end_ORF and end_junction >= start_ORF:
                        condition_array.append(True)
                    else:
                        condition_array.append(False)
                else:
                    condition_array.append(False)
        if sum(condition_array) == 0: outer_condition_array.append(False)
        else: outer_condition_array.append(True)
        outer_condition_detail.append(condition_array)
    
    df['involvement'] = outer_condition_detail
    df['getTranslated'] = outer_condition_array
    df_retained = df[df['getTranslated']==True]
    df_retained = df_retained.drop(columns=['getTranslated']) 
    df_retained = df_retained.set_index(pd.Index(np.arange(df_retained.shape[0])))
    df_filtered = df[df['getTranslated']==False] 
    df_filtered = df_filtered.drop(columns=['getTranslated'])
    df_filtered = df_filtered.set_index(pd.Index(np.arange(df_filtered.shape[0])))
     # filter out events that can not match with existing ones, or events that could match but splicing site won't involve in ORF formation
    return df_retained,df_filtered
            
def alignment_to_uniprot(df,dict_uni_fa,Ens2ACC,mode):
    col1 = []
    col2 = []
    col3 = []
    for i in range(df.shape[0]):
        # collect all uniprot-curated isoform protein sequence
        target = {}
        EnsID = list(df['UID'])[i].split('|')[0].split(':')[1]
        ACCID = Ens2ACC[EnsID] 
        isoforms = list(dict_uni_fa.keys())  # ['Q9NR97','Q9NR97-2'...]
        for iso in isoforms:        
            if ACCID in iso: 
                seq = dict_uni_fa[iso]
                target[iso] = seq
        # collect all mine-predicted isoform protein sequence
        involve = df['involvement'].tolist()[i]  # [True,False,False,True] indicate which transcript would be a good representative
        match_aa = df['ORFaa'].tolist()[i]
        repre = []
        for idx,aa in enumerate(match_aa):
            if involve[idx] == True:   # only consider aa that is caused by splicing event
                bucket = chop_sequence(aa,10)   # chopping
                subnotes = {}
                for j in range(len(bucket)):   # for each 10mer
                    frag = bucket[j]
                    for key,value in target.items():   # for each curated isoform
                        try: 
                            subnotes[key].append(True) if frag in value else subnotes[key].append(False)
                        except KeyError:
                            subnotes[key] = []
                            subnotes[key].append(True) if frag in value else subnotes[key].append(False)
                for k,m in subnotes.items():   # k is key, m is value
                    if sum(m)==0: subnotes[k].append('notAligned')
                    elif sum(m)==len(m): subnotes[k].append('totallyAligned')
                    else: subnotes[k].append('partiallyAligned')
                repre.append(subnotes)
            elif involve[idx] == False:
                repre.append('Either splicing site is not involved in ORF formation or splicing event doesn\'t occur in this certain transcript')
        col1.append(repre)
        # define what kind of peptide this splicing sites would generate by interogratting each repre list
        identity = []
        for n in repre:
            definition = ''
            if isinstance(n,dict):
                for p in n.values():   #n will be {'P14061':[True,True,False,'partiallyAligned'],'P14061-2':[True,True,False,'partiallyAligned']}
                    if p[-1] == 'totallyAligned':  
                        definition = 'one of already documented isoforms'
                        break
            else: 
                definition = 'Either splicing site is not involved in ORF formation or splicing event doesn\'t occur in this certain transcript'
            if definition == '': definition = 'novel isoform'
            identity.append(definition)
        col2.append(identity)
        # let's see if it is possible to generate any novel isoform


        if mode=='strigent':   # as long as one of possible concatenation results in totallyAligned, then we stop pursuing

            for idx,w in enumerate(identity):
                crystal = True   
                if w == 'one of already documented isoforms': 
                    crystal = False
                    break
            if crystal == True:
                crystal_ = []
                for idx,w in enumerate(identity):
                    if w == 'novel isoform': 
                        query_aa = match_aa[idx]
                        try: result = TMHMM(query_aa,'{0}_{1}_{2}'.format(i,idx,EnsID))  # there's some rare case, it can match up with exising one,
                        # but it doesn't have reasonable ORFs, so the returned ORF prediction is ''
                        except: result = False
                        if result:
                            crys = (True,idx)  # we need look into that
                            crystal_.append(crys)
                if crystal_ == []: col3.append(False) # no need to look into this event anymore
                else: col3.append(crystal_)
            else:
                col3.append(False)
            
        elif mode=='loose':    # consider every possible novel isoform possibilities
            crystal = []
            for idx,w in enumerate(identity):
                if w == 'novel isoform': 
                    query_aa = match_aa[idx]
                    try: result = TMHMM(query_aa,'{0}_{1}_{2}'.format(i,idx,EnsID))  # there's some rare case, it can match up with exising one,
                    # but it doesn't have reasonable ORFs, so the returned ORF prediction is ''
                    except: result = False
                    if result:
                        crys = (True,idx)  # we need look into that
                        crystal.append(crys)
            if crystal == []: col3.append(False) # no need to look into this event anymore
            else: col3.append(crystal)
            
        
    df['alignment'] = col1
    df['identity'] = col2
    try:df['interest'] = col3
    except: 
        print(col3)
        raise Exception('hi')
    return df

    
          

    
def read_uniprot_seq(path):
    dict_fa = {}
    with open(path,'r') as in_handle:
        for title,seq in SimpleFastaParser(in_handle):
            uniID = title.split('|')[1]
            dict_fa[uniID] = seq
    return dict_fa       

def chop_sequence(seq,kmer):   # how to splice sequence, elegant way to use range
    frag_bucket = []
    for i in range(0,len(seq),kmer):
        try:
            frag_bucket.append(seq[i:i+kmer])
        except:
            frag_bucket.append(seq[i:])
    return frag_bucket


    
def IDmappingACC2Ensembl(lis,mannual=False):  
    if os.path.exists(os.path.join(outFolder,'Ens2ACC.p')): 
        with open(os.path.join(outFolder,'Ens2ACC.p'),'rb') as f2:
            Ens2ACC = pickle.load(f2)
    else:    
        query = ' '.join(lis)
        url = 'https://www.uniprot.org/uploadlists/'    
        params = {
        'from': 'ACC+ID',
        'to': 'ENSEMBL_ID',
        'format': 'tab',
        'query': query
        }
        
        data = urllib.parse.urlencode(params)  # convert dict to string
        data = data.encode('utf-8')     # convert string to byte
        req = urllib.request.Request(url, data)
        with urllib.request.urlopen(req) as f:
           response = f.read()   
        a = response.decode('utf-8')        # convert byte to string
        
        rows = a.split('\n')
        rows = rows[1:]   # the first item is 'from':'to'
        rows = rows[:-1]  # the last item is empty, it is becasue the last '\n' got recognized as newline
        ACC2Ens,Ens2ACC = {},{}
        for item in rows:
            ACC = item.split('\t')[0]
            Ens = item.split('\t')[1]        
            ACC2Ens[ACC] = Ens
            Ens2ACC[Ens] = ACC
        diff = len(set(lis)) - len(set(list(ACC2Ens.keys())))
        if mannual == False:
            print('There are {0} membrane protein missing!!!!'.format(diff))
        elif mannual == True:
            print('There are {0} membrane protein needs to be mannually checked, you can exit anytime'.format(diff))            
            for acc in lis:     # some ACC doesn't map to a EnsGID, double check see if every one in human membrane list has been mapped
                try: ACC2Ens[acc]
                except KeyError:
                    print('Please mannually check if {0} has corresponding EnsGID:'.format(acc))
                    cond = input('Does {0} have corresponding EnsGID? (y|n|e):'.format(acc))
                    if cond == 'y':
                        EnsGID = input('Please enter the corresponding EnsGID:')
                        ACC2Ens[acc] = EnsGID
                        Ens2ACC[EnsGID] = acc
                        print('Great!, {0} has been mannually added to dictionary.'.format(EnsGID))
                    if cond == 'n':
                        print('Caveat: {0} is a humen membrane protein, it will not be considered in following analysis'.format(acc))
                    if cond == 'e':
                        break
           
            with open(os.path.join(outFolder,'Ens2ACC.p'),'wb') as f1:
                pickle.dump(Ens2ACC,f1)
    return Ens2ACC    # as complete as possible

def TMHMM(aa,name):
    # download TMHMM linux version, untar it.
    # change shabang of tmhmm and tmhmmformat.pl to the path of perl 5+ you loaded
    # export the path of tmhmm to $PATH, finsh configuration

    # in my use case, save those for convenience
    # perl: /usr/local/perl/5.20.1/bin/perl
    # tmhmm: /data/salomonis2/LabFiles/Frank-Li/python3/TMHMM/tmhmm-2.0c/bin
    if not os.path.exists(os.path.join(outFolder,'TMHMM_temp')): os.makedirs(os.path.join(outFolder,'TMHMM_temp'))
    with open(os.path.join(outFolder,'TMHMM_temp','{}.fasta'.format(name)),'w') as f1:
        f1.write('>peptide_{}\n'.format(name))
        f1.write(aa)
    with open(os.path.join(outFolder,'TMHMM_temp','{}.out'.format(name)),'w') as f2:
        subprocess.run(['tmhmm',os.path.join(outFolder,'TMHMM_temp','{}.fasta'.format(name))],stdout = f2)
    with open(os.path.join(outFolder,'TMHMM_temp','{}.out'.format(name)),'r') as f3:
        next(f3)
        punchline = f3.readline().rstrip('\n').split(' ')
        TMn = int(punchline[-1])
    result = True if TMn > 0 else False
    return result

def diffNovelFromNotInvolved(df):
    col = []

    for i in range(df.shape[0]):
        cond = True
        exam_match_whole_tran = df['exam_match_whole_tran'].iloc[i]   # will be a list already, 1*1
        for item in exam_match_whole_tran:
            if item: 
                cond = False
                break
        col.append(cond)
    df['cond'] = col
    try:df_novel = df[df['cond']]
    except:
        print(len(col),df.shape[0],col)
        df.to_csv(os.path.join(outFolder,'fjdfjd.txt'),sep='\t',index=None)
        raise Exception('jkjk')
    return df_novel

def diffNovelFromNotInvolved_new(df):   # we want novel one but not the one that fail the involvement test just because their junction doesn't contribute to ORF formation,those will be ['skipped'] in the third round column since they are normal splicing event 
    col = []
    for i in range(df.shape[0]):
        cond = df['third_round'].iloc[i]
        if cond==['intron'] or cond==['unrecoverable']: 
            crystal = True
            col.append(crystal)
        elif cond == ['skipped']: 
            crystal = False
            col.append(crystal)
        else: 
            hits = sum([True if j else False for j in cond])
            if hits == 0: crystal = True   # in third round, no match, ['','',''], means trans-splicing but no match or novel ordinal, no match
            else: crystal = False           
            col.append(crystal)
    df['check'] = col
    df_novel = df[df['check']]
    df_novel = df_novel.set_index(pd.Index(np.arange(df_novel.shape[0])))
    return df_novel

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


def convertExonList_pep(df):
    dictExonList = {}
    for i in range(df.shape[0]):
        EnsGID = df.iat[i,0]
        EnsPID = df.iat[i,2]
        exonList = df.iat[i,3]
        try: dictExonList[EnsGID][EnsPID] = exonList
        except KeyError: dictExonList[EnsGID] = {EnsPID:exonList}
        # {EnsGID:{EnsPID:exonList,EnsPID:exonList}}
    return dictExonList

def biotype(df):
    dic = {}
    for i in range(df.shape[0]):
        EnsGID = df.iat[i,0]
        EnsPID = df.iat[i,1]
        Anno = df.iat[i,2]
        try:
            dic[EnsGID][EnsPID] = Anno
        except KeyError:
            dic[EnsGID] = {EnsPID:Anno}
    return dic
    # {EnsGID:{EnsPID:Anno,EnsPID:Anno}}

def getORFaa(df):
    col = []
    for i in range(df.shape[0]):
        ORF = df.iloc[i]['ORF']
        if ORF == ['None']:
            col.append('None')
        else:
            tempArray = []
            for transcript in ORF:
                if not transcript: tempArray.append('')
                else:
                    maxAA = str(Seq(transcript,generic_dna).translate(to_stop=False))
                    tempArray.append(maxAA)
            col.append(tempArray)
    df['ORFaa'] = col
    return df


def list_check(lis):
    if lis==['skipped'] or lis==['intron'] or lis ==['unrecoverable']: cond = False
    else:
        hits = sum([True if item else False for item in lis])
        if hits > 0: cond = True
        elif hits == 0: cond = False
    return cond

def build_sorted_exons(EnsGID,exonlists):  # E1.2|E1.3|E2.3|E3.4
    series = []   # store sorted position for each exon
    start_exon = exonlists.split('|')[0]
    strand = dict_exonCoords[EnsGID][start_exon][1]
    if strand == '+':
        start = dict_exonCoords[EnsGID][start_exon][2]
    else:
        start = dict_exonCoords[EnsGID][start_exon][3]  # negative strand, the right most position will be the start, also the largest number

    # final structure {'EnsID':{E1:[chr,strand,start,end],E2:[chr,strand,start,end]}}   
    exonlist = exonlists.split('|')
    dict_judge = {}
    for j in range(len(exonlist)):
        coords = dict_exonCoords[EnsGID][exonlist[j]]
        strand = coords[1]
        judge = check_exonlist_general(exonlist,j,strand)
        dict_judge[exonlist[j]] = judge
    
        
    dic = {}
    for subexon in exonlist:
        exon_num = subexon.split('.')[0]
        subexon_num = subexon.split('.')[1]
        if exon_num in dic:
            dic[exon_num].append(subexon_num)
        else:
            dic[exon_num] = []
            dic[exon_num].append(subexon_num)  # E14 > [1,2,4,5]
    accum = 0
    for exon,sub in dic.items():
        incre,position = check_consecutive(exon,sub,dict_judge,EnsGID,strand,accum)
        accum += incre
        series.extend(position)
    series.sort()   # ascending order [5,9,15,...]
    series = [0]+series
    return series
        
            
                    

            


def check_consecutive(exon,sub,dict_judge,EnsGID,strand,accum):   # E14  > [1,2,4,5]
    #print(exon,sub,dict_judge,accum)
    position = []
    lis_int = [int(x) for x in sub]
    diff1 = np.diff(lis_int,1)   # array([1,2,1])
    diff1 = [int(x)-1 for x in diff1]    # [0,1,0]
    split = np.nonzero(diff1)[0].tolist()  # if pos=1, it means in original list, after index 1 will have a breaking point
    #print(split)
    if split:   # have breaking point
        split = [y + 1 for y in split]    
        # lis_int contains original list, split contains all the indices that identical to the first one in each subgroup
        result=[lis_int[i:j] for i,j in zip([0]+split,split+[None])] 
        for chunk in result:  # chunk[1,2], chunk[4,5]
            query_s = str(exon)+'.'+str(chunk[0])
            query_e = str(exon)+'.'+str(chunk[-1])
            if strand=='+':
                start = dict_exonCoords[EnsGID][query_s][2]
                end = dict_exonCoords[EnsGID][query_e][3] if dict_judge[query_e] else int(dict_exonCoords[EnsGID][query_e][3])-1
                relaPos = int(end) - int(start) + 1   # think 5-1=4, but 5 will be the 5th one
                position.append(relaPos+accum)
            elif strand == '-':
                start = dict_exonCoords[EnsGID][query_s][3]
                end = dict_exonCoords[EnsGID][query_e][2] if dict_judge[query_e] else int(dict_exonCoords[EnsGID][query_e][2])+1
                relaPos = int(start) - int(end) + 1
                position.append(relaPos+accum)
    else:   # E15 > [1,2,3]   3 consecutive 
        query_s = str(exon) + '.' + str(sub[0])
        query_e = str(exon) +'.'+ str(sub[-1])
        #print(query_s,query_e)
        
        if strand=='+':
            start = dict_exonCoords[EnsGID][query_s][2]
            end = dict_exonCoords[EnsGID][query_e][3] if dict_judge[query_e] else int(dict_exonCoords[EnsGID][query_e][3])-1
            relaPos = int(end) - int(start) + 1   # think 5-1=4, but 5 will be the 5th one
            position.append(relaPos+accum)
        elif strand=='-': 
            start = dict_exonCoords[EnsGID][query_s][3]
            end = dict_exonCoords[EnsGID][query_e][2] if dict_judge[query_e] else int(dict_exonCoords[EnsGID][query_e][2])+1
            relaPos = int(start) - int(end) + 1
            position.append(relaPos+accum)
        #print(relaPos)

                
    return relaPos,position


def check_translation(EnsGID,EnsPID):
    if EnsPID == 'None':    # usually from RefSeq dataset
        result = '*'
    elif 'PEP' in EnsPID:   # ENSP854949-PEP
        result = '*'
    else:   # protein coding gene or NMD
        pepAnno = dict_biotype[EnsGID]  #{ENSP:anno,ENSP:anno}
        if pepAnno[EnsPID] == 'protein_coding': result = '#'
        else: result = '*'
    return result
    
        






# https://rest.ensembl.org/documentation/info/lookup
def grabEnsemblTranscriptTable(EnsTID):

    server = "https://rest.ensembl.org"
    ext = "/lookup/id/{0}?expand=1".format(EnsTID)     
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})     
    try: decoded = r.json()
    except: 
        print('JSON unknoen error')
        result = '#'     # I don't think if running on local, this condition will ever be reached
    else:
        try:
            biotype = decoded['biotype']
        except:
            result = '*'   # unknown crash, might be a invalid Ensembl ID

        if biotype == 'protein_coding': result = '#'
        else: result = '*'   # non-protein coding genes or data from other source
        
    return result

def ORF_check(df):
    col1,col2 = [],[]
    for i in range(df.shape[0]):
        print('The {}th run'.format(i))
        temp=uid(df,i) 
        EnsGID = list(temp.keys())[0].split(':')[1]
        space = dictExonList_p[EnsGID]  # [('ENSP',exonlists),('ENSP',exonlists)...]
        ORF = df.iloc[i]['ORF']
        first = df.iloc[i]['exam_first_whole_transcripts']
        second = df.iloc[i]['second_round']
        third =df.iloc[i]['third_round']
        if list_check(first): whole = first
        elif not list_check(first) and list_check(second): whole = second
        elif not list_check(first) and not list_check(second) and list_check(third): whole = third
        else:
            NMD = ['None']   # intron, unrecoverable or third round is still all empty
            translate = ['None']
            col1.append(NMD)
            col2.append(translate)
            continue

        NMD = []
        translate = []
        #print(len(ORF),len(space))
        if len(ORF) == len(space):   # not trans-splicing events
 
            for j in range(len(ORF)):
                orf = ORF[j]
                if not orf: 
                    NMD.append('')
                    translate.append('')
                elif orf:
                    whole_ = whole[j]
                    try:space_ENSP = space[j][0]
                    except:print(j,space,EnsGID,ORF); print(df);raise Exception
                    
                    space_exons = space[j][1]
                    result = check_translation(EnsGID,space_ENSP)

                    translate.append(result)
                    
                    series = build_sorted_exons(EnsGID,space_exons)
                    num_exon = len(series) - 1
                    #print(series,num_exon)
                    #print(orf,type(orf))
                    orf_end_pos = whole_.find(orf)+len(orf)-1
                    residing = bisect.bisect_left(series,orf_end_pos)   # which exon it resides on
                    #print(residing)
                    if residing <= num_exon-2: NMD.append('*')   # potentially NMD
                    else: NMD.append('#')  # good candidate
        else:   # trans-splicing events
            for j in range(len(ORF)):
                orf = ORF[j]
                if not orf: 
                    NMD.append('')
                    translate.append('')
                elif orf : 
                    if orf=='None': 
                        NMD.append('None')
                        translate.append('None')
                    else:
                        NMD.append('#')    # currently don't support interrogation of NMD for tran-splicing event
                        translate.append('#')   
        
        col1.append(NMD)
        col2.append(translate)
    df['NMD_check'] = col1
    df['translate_check'] = col2
    return df

def getORF(df):
    col = []
    for i in range(df.shape[0]):
        first_round = df.iloc[i]['exam_first_whole_transcripts']
        second_round = df.iloc[i]['second_round']
        third_round = df.iloc[i]['third_round']
        if third_round == ['intron'] or third_round == ['unrecoverable']:
            col.append(['None'])
        elif not third_round == ['skipped']:   # transcripts are in third_round
            tempArray = []
            for transcript in third_round:
                if not transcript: tempArray.append('')
                else:
                    transcript = transcript.replace(',','')
                    maxTran = transcript2peptide(transcript)
                    tempArray.append(maxTran)
            col.append(tempArray)        # so ['','',''] ORF could be either unrecoverable third round or the predicted ORF is too short
        elif third_round == ['skipped']:
            if not second_round == ['skipped']:   # transripts are in second_round
                tempArray = []
                for transcript in second_round:
                    if not transcript: tempArray.append('')
                    else:
                        maxTran = transcript2peptide(transcript)
                        tempArray.append(maxTran)
                col.append(tempArray)  
            elif second_round == ['skipped']:   # transcripts are in first_round
                tempArray = []
                for transcript in first_round:
                    if not transcript: tempArray.append('')
                    else:
                        maxTran = transcript2peptide(transcript)
                        tempArray.append(maxTran)
                col.append(tempArray) 
    df['ORF'] = col
    return df

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
        if re.search(rf'{re.escape(Exons1)}',item) or re.search(rf'{re.escape(Exons2)}',item) or re.search(rf'{re.escape(query)}$',item):
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
                judge_query = dict_judge[exam2]
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
                Remember: the arguments to query_from_dict_fa is very simple, it is just the coord[2] and coord[3],
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


def third_round(df_second):  # after second run
    '''
    In second round
    1. [''] means trans-splicing(non-trailing), novel ordinal, intron retention, they don't have '_' in the exon
    2. ['','','','',...''] means newsplicingsite, trans-splicing(trailing) or Alt5,3 but even trimming the trailing part can not match them to existing ones    
    
    '''
    col = []
    for i in range(df_second.shape[0]):
        second = df_second.iloc[i]['second_round']
        temp=uid(df_second,i) 
        EnsGID_this = list(temp.keys())[0].split(':')[1]
        exam1 = list(temp.values())[0][0].split('-')[0]   # E22
        exam2 = list(temp.values())[0][0].split('-')[1]   # ENSG:E31
        if second == [''] and 'ENSG' in exam2:  # trans-splicing(non_trailing)
            EnsGID_trans = exam2.split(':')[0]
            exam2_trans = exam2.split(':')[1]
            full_left = single_left_match(exam1,EnsGID_this)
            full_right = single_right_match(exam2_trans,EnsGID_trans)
            full = cat_left_right_asym(full_left,full_right)
            col.append(full)
        elif 'I' in (exam1 + exam2) and second == ['']:   # intron retention
            col.append(['intron'])
        elif second == ['']:   # novel ordinal    # E3.4(exam1) - E5.1(exam2)
            full_left = single_left_match(exam1,EnsGID_this)
            full_right = single_right_match(exam2,EnsGID_this)
            full = cat_left_right_sym(full_left,full_right)
            col.append(full)
        else: # ['skipped'] or ['TTTTAAA','AATTGGCC'] has matches or ['','',''] that we don't want to recover
            hits = sum([True if item else False for item in second])
            if hits > 0: col.append(['skipped'])
            elif hits == 0: col.append(['unrecoverable'])  # means point 2 in comments above
    df_second['third_round'] = col  # if third_round is still ['','',''], means trans and novel ordinal, their single subexon still don't match up with any
    return df_second

            
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
            if i and j: result.append(i + ',' + j)
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
    
    
def second_round(df_first):
    col = []
    for i in range(df_first.shape[0]):  
        first = df_first.iloc[i]['exam_first_whole_transcripts']
        hits = sum([True if match else False for match in first])
        if hits > 0: col.append(['skipped'])   # already matched in first round, not a novel event
        else:
            temp=uid(df_first,i)
             #{'gene:ENSid':[E22-E33,E34-E56]}
             # if fusion gene: E22-ENSG:E31
            EnsID=list(temp.keys())[0].split(':')[1]
            exam1 = list(temp.values())[0][0].split('-')[0]
            exam2 = list(temp.values())[0][0].split('-')[1]
            if '_' in exam1 and not '_' in exam2:    # mode 1
                exam1_exon = exam1.split('_')[0]
                exam1_coord = exam1.split('_')[1]
                query = exam1_exon + '|' + exam2
                result = second_match(EnsID,query,exam1_coord=exam1_coord)
                col.append(result)

            if '_' not in exam1 and '_' in exam2:   # mode 2
                exam2_exon = exam2.split('_')[0]
                exam2_coord = exam2.split('_')[1]
                query = exam1 + '|' + exam2_exon
                result = second_match(EnsID,query,exam2_coord=exam2_coord)
                col.append(result)
                
            if '_' in exam1 and '_' in exam2:
                exam1_exon = exam1.split('_')[0]
                exam1_coord = exam1.split('_')[1]                
                exam2_exon = exam2.split('_')[0]
                exam2_coord = exam2.split('_')[1]
                query = exam1_exon + '|' + exam2_exon
                result = second_match(EnsID,query,exam1_coord=exam1_coord,exam2_coord=exam2_coord)
                col.append(result)
            if not '_' in exam1 and not '_' in exam2:
                result = ['']         # novel ordinal and intron retention and tran-splicing(non-trailing)
                col.append(result)
    df_first['second_round'] = col
    return df_first

def extract_EnsID(df):
    
    UID = list(df['UID'])
    EnsID_array = []
    for item in UID:
        EnsID = item.split('|')[0].split(':')[1]
        EnsID_array.append(EnsID)
    return EnsID_array

def matchWithExonlist(df,df_exonlist,dict_exonCoords):
    col1,col2 = [],[]        
    for i in range(df.shape[0]):
        temp=uid(df,i)
        EnsID=list(temp.keys())[0].split(':')[1]
        Exons_examined = exon_extract(temp,0,EnsID)
        Exons_back = exon_extract(temp,1,EnsID)
        col1.append(core_match(df_exonlist,dict_exonCoords,EnsID,Exons_examined))
        col2.append(core_match(df_exonlist,dict_exonCoords,EnsID,Exons_back))
    df['exam_first_whole_transcripts'] = col1
    #df['back_first_whole_transcripts'] = col2
    return df

def main(intFile,dataFolder,outFolder,mode,cutoff_PSI,cutoff_sampleRatio,cutoff_tissueRatio):
#    intFile = parser.intFile
#    dataFolder = parser.dataFolder
#    outFolder = parser.outFolder
#    print(intFile,dataFolder,outFolder)
    # get increased part
    df_ori = pd.read_csv(intFile,sep='\t')
    #df_ori = GetIncreasedPart(df)


    
    # load the files
    global df_exonlist
    global dict_exonCoords
    global dict_fa
    global dictExonList
    global dictExonList_p
    global df_biotype
    global dict_biotype
    global dicTissueExp

    df_exonlist = pd.read_csv(os.path.join(dataFolder,'mRNA-ExonIDs.txt'),sep='\t',header=None,names=['EnsGID','EnsTID','EnsPID','Exons'])
    dict_exonCoords = exonCoords_to_dict(os.path.join(dataFolder,'Hs_Ensembl_exon.txt'),'\t')   
    dict_fa = fasta_to_dict(os.path.join(dataFolder,'Hs_gene-seq-2000_flank.fa'))
    dictExonList = convertExonList(df_exonlist)
    dictExonList_p = convertExonList_pep(df_exonlist)
    df_biotype = pd.read_csv(os.path.join(dataFolder,'Hs_Ensembl_transcript-biotypes.txt'),sep='\t')
    dict_biotype = biotype(df_biotype)
    print('Please wait 10 minutes for loading GTEx data')

    if not os.path.exists(os.path.join(outFolder,'after_GTEx.txt')):

        with bz2.BZ2File(os.path.join(dataFolder,'dicTissueExp.pbz2'),'rb') as f1:
            dicTissueExp = cpickle.load(f1)  
        print('Already loaded GTEx data')


        # get tumor-specific part
        df_ori =  check_GTEx(df_ori,cutoff_PSI,cutoff_sampleRatio,cutoff_tissueRatio)   
        if df_ori.shape[0] == 0: 
            raise Exception('After checking GTEx, no events remain')
        df_ori.to_csv(os.path.join(outFolder,'after_GTEx.txt'),sep='\t',index=None)
    
    else:
        df_ori = pd.read_csv(os.path.join(outFolder,'after_GTEx.txt'),sep='\t')

    print('Overlapping with human membrane proteins')
    
    # overlapping with human membrane proteins
    df_membraneProteins = pd.read_csv(os.path.join(dataFolder,'human_membrane_proteins.txt'),sep='\t')  
    ACClist = df_membraneProteins['Entry'].tolist()
    Ens2ACC = IDmappingACC2Ensembl(ACClist,False) 
    EnsID = extract_EnsID(df_ori)
    df_ori['condition'] = [True if item in list(Ens2ACC.keys()) else False for item in EnsID]
    df_ori_narrow = df_ori[df_ori['condition'] == True]
    if df_ori.shape[0] == 0:
        raise Exception('After overlapping with membrane proteins, no events remain')
    df_ori_narrow = df_ori_narrow.drop(columns=['condition'])
    df_ori_narrow = df_ori_narrow.set_index(pd.Index(np.arange(df_ori_narrow.shape[0])))



    # match with all existing transcripts
    print('first round matching...')
    df_first = matchWithExonlist(df_ori_narrow,df_exonlist,dict_exonCoords)  
    print('second round matching...')
    df_second = second_round(df_first)
    print('third round matching...')
    df_third = third_round(df_second)
    print('predicting most likely ORFs...')
    df_ORF = getORF(df_third)
    #df_ORF.to_csv(os.path.join(outFolder,'inter.txt'),sep='\t')
    print('labelling potential ORFs that will subjected to NMD and non-translatable ones...')
    df_ORF_check = ORF_check(df_ORF)
    print('in-silico translation...')
    df_ORF_aa = getORFaa(df_ORF_check)  
    
    # derive the junction site sequence and add two columns to df_ori
    new_df_narrow = retrieveJunctionSite(df_ORF_aa,dict_exonCoords,dict_fa)
    
    df_retained,df_filtered = check_if_good_representative(new_df_narrow)
    
    # for retained ones

    
    ### final alignment
    dict_uni_fa = read_uniprot_seq(os.path.join(dataFolder,'uniprot_isoform.fasta'))
        
    df_retained_aligned = alignment_to_uniprot(df_retained,dict_uni_fa,Ens2ACC,mode) 
    df_retained_aligned = df_retained_aligned.drop(columns=['exam_seq'])  
    df_retained_aligned.to_csv(os.path.join(outFolder,'df_retained.txt'),sep='\t',index=None)
    
    # for filtered ones
    df_novel = diffNovelFromNotInvolved_new(df_filtered)
    #df_novel = df_novel.drop(columns=['exam_seq'])
    df_novel.to_csv(os.path.join(outFolder,'df_novel.txt'),sep='\t',index=None)    

def usage():
    print('Usage:')
    print('python3 NeoEpitopePredictor.py --intFile /data/salomonis2/LabFiles/Frank-Li/receptor/breast/pseudoPSImatrix.txt --dataFolder /data/salomonis2/LabFiles/Frank-Li/python3/data --outFolder /data/salomonis2/LabFiles/Frank-Li/receptor/breast --mode strigent --cutoffPSI 0.1 --cutoffSample 0.1 --cutoffTissue 0.1')
    print('Options:')
    print('--intFile : path of input file')
    print('--dataFolder : path of data folder')
    print('--outFolder : output folder')    
    print('--mode : using TMHMM or not')
    print('--cutoffPSI : above this PSI value will be labelled as expressed')
    print('--cutoffSample : above this ratio this tissue will be labelled as expressed')
    print('--cutoffTissue : above this ratio this event will be labelled as expressed in normal tissue in general')
    
    

if __name__ == "__main__":
    #os.chdir('/Users/ligk2e/Desktop/project_LUAD')
    
    try:
        options, remainder = getopt.getopt(sys.argv[1:],'h',['help','intFile=','dataFolder=','outFolder=','mode=','cutoffPSI=','cutoffSample=','cutoffTissue='])
    except getopt.GetoptError as err:
        print('ERROR:', err)
        usage()
        sys.exit(1)
    for opt, arg in options:
        if opt in ('--intFile'):
            intFile = arg
            print('Input file is:', arg)
        elif opt in ('--dataFolder'):
            dataFolder = arg
            print('Data folder is:',arg)
        elif opt in ('--outFolder'):
            outFolder = arg
            print('output folder:',arg)
        elif opt in ('--mode'):
            mode = arg

            print('Using TMHMM?:',arg)
        elif opt in ('--cutoffPSI'):
            cutoff_PSI = float(arg)
            print('cutoff for PSI value:',arg)
        elif opt in ('--cutoffSample'):
            cutoff_sample = float(arg)
            print('cutoff for sample ratio:',arg)    
        elif opt in ('--cutoffTissue'):
            cutoff_tissue = float(arg)
            print('cutoff for tissue ratio:',arg)        
        elif opt in ('-h','--help'):
            usage()
            sys.exit(1)

    main(intFile,dataFolder,outFolder,mode,cutoff_PSI,cutoff_sample,cutoff_tissue)
    
    
    
    
#    parser = argparse.ArgumentParser(description='Receptor Protein Prediction') # ArgumentParser object
#    parser.add_argument('--intFile',type=str,default='.',help='input file path')
#    parser.add_argument('--dataFolder',type=str,default='./data',help='data folder path')
#    parser.add_argument('--outFolder',type=str,default='.',help='output folder path')
#    args = parser.parse_args()   # namespace object
#    main(args)


    
    # summarize the distribution of splicing event in df_all and df_increased
#    freq_all = ChroDistribution(df)
#    freq_increased = ChroDistribution(df_ori)
#    print(freq_all,freq_increased)
    
    # continue exploit on chromosomes
#    chro_dict = {
#            'chr1': [1961,248,'Metacentric'],    #[1961genes,248 or so million bp, type of centromere]
#            'chr2': [1194,242,'Submetacentric'],
#            'chr3': [1024,198,'Metacentric'],
#            'chr4': [727,190,'Submetacentric'],
#            'chr5': [839,181,'Submetacentric'],
#            'chr6': [996,170,'Submetacentric'],
#            'chr7': [862,159,'Submetacentric'],
#            'chr8': [646,145,'Submetacentric'],
#            'chr9': [739,138,'Submetacentric'],
#            'chr10': [706,133,'Submetacentric'],
#            'chr11': [1224,135,'Submetacentric'],
#            'chr12': [988,133,'Submetacentric'],
#            'chr13': [308,114,'Acrocentric'],
#            'chr14': [583,107,'Acrocentric'],
#            'chr15': [561,101,'Acrocentric'],
#            'chr16': [795,90,'Metacentric'],
#            'chr17': [1124,83,'Submetacentric'],
#            'chr18': [261,80,'Submetacentric'],
#            'chr19': [1357,58,'Metacentric'],
#            'chr20': [516,64,'Metacentric'],
#            'chr21': [215,46,'Acrocentric'],
#            'chr22': [417,50,'Acrocentric'],
#            'chrX': [804,156,'Submetacentric'],
#            'chrY': [63,57,'Acrocentric']}  
#    PlotChroScarse(chro_dict,'human chromosome genes distribution.pdf')

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    