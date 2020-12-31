# AltNeo-T Full Documentation (Users and Developers)

AltNeo-T is an automated program to prioritize **T**umor **S**pecific **A**ntigen (**TSA**), aka **Neoantigen**, generated from abnormal alternative splicing events, which can be presented by MHC molecule on cell surface and trigger T cell response. Hence, we can utilize this tool to identify potential cancer immunotherapy target by which we can employ to boost patients' own cellular immunity. Please check its pair program AltNeo-B for B cell immunity.

AltNeo-T takes the input of a list of **L**ocal **S**plicing **E**vents(**LSV**). Theoritically it can accept LSV list from any alternative splicing detection program like MultiPath-PSI, rMATs, MAJIQ, LeafCutter, etc. Currently we develop the pipeline specifically for **MultiPath-PSI**. For detail of MultiPath-PSI, please refer to its documentation page (https://altanalyze.readthedocs.io/en/latest/Algorithms/#multipath-psi-splicing-algorithm). In the next section, I will illustrate how this AltNeo-T program work.


## Input File

Let's begin with the input file, it would be something like this:

UID | Optional
--- | --- 
LHPP:ENSG00000107902:E10.1-E12.1 | Ordinal
UNC93B1:ENSG00000110057:E5.1-E6.2_67996641 | Alt-3
NFIX:ENSG00000008441:I40.1_13076665-E41.1 | Novel Exon
SF3B1:ENSG00000115524:I4.1-E5.1 | Intron Retention
HBG2:ENSG00000196565:E14.2-ENSG00000213934:E3.1 | Trans-splicing
RNF123:ENSG00000164068:U0.1_49689185-E2.1 | UTR event
SYNGR1:ENSG00000100321:E7.1_39364266-E8.1 | Alt-5

When you actually run the program, the input LSV would be the one shown above in conjunction with a "background event (LSV)". For brevity, let's stick the the aforementioned notation for now. The input file can be any **tab delimited format with a minimum requirement that there must exist one column named "UID"**, you can specifiy as many additional columns as you want, which are totally optional.

**Figure 1**

![Figure1: whole workflow](https://github.com/frankligy/AltNeo-BT/blob/main/images/Figure1.png)

## Loading necessary files 

AltNeo-T will load a bunch of necesary files from `data` folder you specified via command line parameter. It will contains:

1. **mRNA-exonID.txt** (Each ENST ID and its subexons constitutions)
2. **Hs_Ensembl_exon.txt** (Each subexon and its associated attributes, namely, chromosome, strand, start_pos, end_pos, detail shown in **Figure 2**) There is an issue with the subexon coordinates called **overlapping issue**, it is illustrated in **Figure 3**
3. **Hs_gene-seq-2000_flank.fa** (Sequences of each ENSG, if the gene lies on reverse strand, the sequences stored in this file is reverse strand sequence, also, there is 2000 bp **buffer sequence** in front of eachc sequence, it means, if the coordinates of a certain ENSG is shown as 50|100, then the first necleotide of this gene lies on the 2001st letter in the stored sequence)

4. **gtfEnsembl91.txt** (the start exon coordinates for each ENST, will be useful when mannually derive intron events translational phase and Nmer from their junction sequences)

**Figure 2**

![Figure2: understanding subexon coordinate](https://github.com/frankligy/AltNeo-BT/blob/main/images/Figure2.png)

**Figure 3**

![Figure3: overlapping issues](https://github.com/frankligy/AltNeo-BT/blob/main/images/Figure3.png)

## GTEx check

We want to focus on neoantigen so it is vital to remove any splicing events that express in normal tissue in an appreciable level, their resultant neoantigens can not be used in clinical therapy because it will harm the normal cells. To this end, we investigate each LSV's expression level in GTEx dataset. The GTEx raw fastq files have been pre-processed to two files (`GTEx_EventAnnotation.txt` and `GTEx_SRARunTable.txt`) and I've generated a hdf5 file (`dicTissueExp.hdf5`) and a pbz2 comparessed file (`dicTissueExp.pbz2`) which store a huge two-layer dictionary whose first level key splicing event, second level key is tissue and then the expression for each event in each tissue at each sample was stored as an array. The code for generating these two files are in `queryGTEx_PSI.py` If you want to query single event, hdf5 file woule be faster since it doesn't need to load whole dictionary into RAM, instead, if you want to query tens of thousands event, please load pbz2 file into RAM first. This step will output an intermediate files which contains Neojunctions that pass the GTEx check. The downstream analysis will perform upon these Neojunctions.

Before GTEx check, we instantiate `Meta` object, then after GTEx check, we instantiate `NeoJ` object, we are ready to spawn subprocess using python `multiprocessing` package `Pool` function. Basically, it divides the huge input dataframe/series into `n_cores` chunks, each chunk will have exactly same structure as the original ones but less rows. For each chunk, we will start our processing.

We also developed a visulization html page called GTEx Viewer which can generate all expression information for a certain queried splicing event. Just simply run `queryGTEx_PSI.py`

We will also be able to report a continuous specificity score to quantify how tumor-specific each splicing event is.

## Retrieve Junction sequence

As the first step, we retrieve the junction sequence from the splicing event. The logic is straightforward, you have the subexon so you can obtain their chromosomal coordinates, you use these coordinates to retrieve sequence from fasta files. But different types of splicing evnets would need differnt care, the most important function for this step is `subexon_tran`, which returns the sequence of a queried subexon, no matter what type it is. `subexon_tran` will call a low-level function `query_from_dict_fa` to actually retrieve the sequence from fasta file, here the arguments for `query_from_dict_fa` is very simple, **it is just the forward strand coordinates of the start and end position, for reverse strand subexon, it is the forward strand coordinates of its end and start position.**

If it is a UTR event, in `subexon_tran` function, it will call for `utr_junction` function which will establish a GET request to USCS genome browser to retrieve the UTR sequence, we only cut UTR sequence of length 100. **This step needs internet connection**. 

## Match each LSV to all potential parantal transcripts, and get whole transcripts' sequence

We want to engraft each LSV event to all of its originated parental transcripts. You may wonder why? It is because if we only use junction sequence and do in-silico translation in a 3-way fashion, it will generate a lot of false positive neoantigen because in reality, only one of 3 translational frame will be adopted. Then how can we infer this correct translational frame? We want to use the whole transcript information, if we know the **T**ranslational **S**tart **S**ite (**TSS**) and the region that junction sequence lies on the whole transcript, we can accurately derive the translational frame, or we call phase. The phase can be in the value of 0, 1 or 2. The way I defined them is shown in **Figure 4**


**Figure 4**

![Figure4: translational frame](https://github.com/frankligy/AltNeo-BT/blob/main/images/Figure4.png)

Now we know the reason why we bother to match LSV to its pareantal transcript. The way we actually implement it is to employ a three round search. In the first round search, carried out by function `matchWithExonlist` and its child function `core_match`, we match all the ordinal event to its parental transcript. The idea is fairly simple, LSV `E1.2|E3.1` will match to transcript `E1.1|E1.2|E3.1|E4.1` right? That's it. We just map each LSV event to all its documented transcripts, see if it can match. Hence, the result will be a list of length equal to all documented transcripts for a certain gene (ENSG), if it successfully match the transcript, we will return the whole transcript cDNA sequence, if not, it will return a empty string "". So when you inspect the result from this step, you will see ["","AAAAC...CCTT","TTTT...CCC"], now we know what it means. A rare situation that could happen is that it return a empty list [], it means the ENSG we queried doesn't exist in `mRNA-exonID.txt` file. 

While we've already mapped all ordinal event to its parental transcripts, what about those empty string ""? They may denote that a LSV `E1.2|E3.1` really doesn't fit in a transcript `E1.1|E3.2|E4.1`, but the another chances are the LSV looks like this `E1.2_483943993|E3.1` which can not map to trancript `E1.1|E1.2|E3.1|E4.1` directly, but it is highly likely that this LSV is derived from this trancript. You initiate a second round search to recover all the Alt-5 and Alt-3 type events. It is carried out by functions `rescueEvent_1` `second_round_mhc` `second_match`. The idea behind this step is also very easy to understand, let's explain it in the next paragraph.

We first categorize the splicing events we want to rescue into three modes: 

1. mode1: `E1.2_483943993|E3.1`
2. mode2: `E1.2|E3.1_4374343438`
3. mode3: `E1.2_483943993|E3.1_4374343438`

No matter which mode it is, we will remove the trailing part and all three will degrate to `E1.2|E3.1`, we use these degenerated version to map to all transcripts. Now it is able to match up with `E1.1|E1.2|E3.1|E4.1`. We put `E1.1` as variable `bucket_left`, `E3.1|E4.1` as variable `bucket_right`, and we will retrieve the genomic sequence for `bucket_left` and `bucket_right` individually. Then the most critical part is the `E1.2`, we will directly get its sequence `query_frag` by calling `query_from_dict_fa` function. This part is what these three modes will differ from each othter, because of the overlapping issue I meantioned before, the arguments we pass to function `query_from_dict_fa` will differ depending on the mode and the strand, it can be easily illustrated through **Figure 2** and corresponding codes. Finally, we just simply concatenate `bucket_left`, `query_frag`, `bucket_right` together as the final resultant sequence.

(Notes for me only, the query_frag for mode2 and mode3 may have bugs, reverse strand situtation, we should use start coordinates instead of end coordinates)

We only initiate second round when the resultant list from round one contains zero sequence but all empty strings. After second round, some will return [""], it means there's no "_" in neither exon1 (exon before the vertical line "|") nor exon2 (exon after the vertical line "|"). It signifies (1) trans-splicing event without trailing part, (2) novel ordinal event (event that is ordinal but can not map to any existing transcripts), (3) Intron retention. 

It may also return result like this ["","",""], it means the events that (1) trans-splicing with trailing part, or (2) Alt-5 or Alt-3 events that can not map to existing ones even we trim off the trailing part. 

In the third round, we aim to rescue two types of events, (1) trans-splicing without trailing part, (2) novel ordinal event, also we will label Intron retention event as "intron" and the remain ones as "unrecoverable". The logics behind third round is like this, for trans-splicing event, we individually map exon1 and exon2 (definition is in previous paragraph), and then use `cat_left_right_asym` to put them together, it is asynchronous because we don't require left return and right return list have same length since they come from two different genes which have different number of existing transcripts. For novel ordinal events, we take the same strategy but use `cat_left_right_sym` function to do concatenation, where we require that one have to have matched part for both left and right to allow the concatenation to occur.

It is only a note to myself, I think in this way the novel exon will be recovered by second round but it will generate the whole intron region, which is not correct.

## ORF Finding

For each attained whole trancript sequence, we want to find the optimal or most likely ORF from that, because eventually what we care about the TSS. We developed an customized algorithm to achieve this goal. Frist we use human canonical start codon and stop codon as anchor to cut a bunch of possible ORFs. Then we rank them based on three criteria:

1. Length
2. GC content
3. Coding bias/potential

Length will have top priority under the assumption that longer sequence hold higher chance to translate to genuine peptides. But if two ORFs' length only differ in the range of [0,8], we will treat them as same length and further decide which one is optimal by comparing their GC content and coding potential. GC content is computed as the total occurence of G and C normalized by ORF length, coding potential is the total coding potential (nornalized coding potential by minmaxScalar) for all neighbor triplets normalized by the numebr of triplets in each ORF.

Finally, we will get one optimal ORF for each whole trancript, it might be an empty string if no optimal ORF returned. We will also generate translated peptide from this ORF, but this step is not critical.

## Translation Phase determination

Finally, we reach the goal why we bother deriving the parental trancripts for each LSV. Now we will map the junction sequence, TSS to whole trancript and derive translational phase, illustrated in **Figure 4**. It is noteworthy that to further control the false discoveries, we restrict the splicing events in a way that if the junction region falls out of predicted optimal ORF, we will not consider this splicing event-ORF pair, illustrated in **Figure 5**. If in the previous step, we fail to recover its whole transcript, we may need to trigger "mannual" function to treat them individually and derive tranlational phase, it may be not able to generate perfect phase but we will try our best.

When having translational phase, we can derive Nmer, N can be 8,9,10 and 11. This is the whole search space that will be subjected to MHC presentation prediction (binding prediction).

**Figure 5**

![Figure 5](https://github.com/frankligy/AltNeo-BT/blob/main/images/Figure5.png)


## MHC binding affnity prediction

Since it's a well-established field, we don't want to re-invent wheels, two state-of-art methods, namely netMHCpan-4.1 which is a wrapped software package, and MHCflurry-2.0, which is a python package, were used in AltNeo-T program. It supports both methods. This is the most time-consuming step and may take several hours to finish up. 

## MHC-I immunogenicity prediction

Binding doesn't directly infer immunogenicty, we developed DeepImmuno, a Convolutional Neural Network based approach to predict each neoantien's immunogenicity. Please check out its github repository (https://github.com/frankligy/DeepImmuno) and our preprint (https://www.biorxiv.org/content/10.1101/2020.12.24.424262v1).



## Usage

Each time you can only query neoantigen of a single length, i.e. 9mer, we usually want to have all 8-11mer in MHCI case. We wrote a wrapper shell script `sub_mhcPresent,sh` to wrap the main python script with two other auxiliary python files `get_combined.py` and `get_Augment.py`. Which can expedite this process.

## Contact

This tool is develped by Guangyuan (Frank) Li, a bioinformatics PhD student in Cincinnati Children's Hospital Medical Center (CCMHC) and University of Cincinnati (UC). Don't hesitate to ask questions (li2g2@mail.uc.edu)

