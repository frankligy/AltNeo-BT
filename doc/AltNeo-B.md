## AltNeo-B (For Users and Developers)

AltNeo-B is an automated program for identifying surface proteins in tumor cell whose extracellular domain is different from its counterpart in normal cell, caused by dysregulated alternative splicing events. These aberrant extracellular regions can serve as cancer immunotherapy target for B-cell based therapy, for example, CAR-T therapy. 

It is paired with AltNeo-T program, which aims to identify cancer immunotherapy target for T-cell based therapy, the documentation for AltNeo-T is (https://github.com/frankligy/AltNeo-BT/blob/main/doc/AltNeo-T.md). 

It shares a lot of similarities with AltNeo-T in terms of their input file type and pre-processing, but I will re-iterate them again in this documentation. The input file looks like this, it is just a list of Local Splicing Event (**LSV**) from MultiPath-PSI program.

UID | Optional
--- | --- 
LHPP:ENSG00000107902:E10.1-E12.1 | Ordinary
UNC93B1:ENSG00000110057:E5.1-E6.2_67996641 | Alt-3
NFIX:ENSG00000008441:I40.1_13076665-E41.1 | Novel Exon
SF3B1:ENSG00000115524:I4.1-E5.1 | Intron Retention
HBG2:ENSG00000196565:E14.2-ENSG00000213934:E3.1 | Trans-splicing
RNF123:ENSG00000164068:U0.1_49689185-E2.1 | UTR event
SYNGR1:ENSG00000100321:E7.1_39364266-E8.1 | Alt-5

Now let's begin the journey!

**Figure 1**

![Figure1](https://github.com/frankligy/AltNeo-BT/blob/main/images/B-Figure1.png)

## Load necessary Dataset

Likewise, AltNeo-B will load some necessary datasets for downstream analysis, they are:

1. **mRNA-exonID.txt**
2. **Hs_Ensembl_exon.txt**
3. **Hs_gene-seq-2000_flank.fa**

Please check AltNeo-T documentation for detailed explanation for above three datasets.

4. **Hs_Ensembl_transcript-biotypes.txt** (This file stores the biotype annotation for each ENSP, it will be used when doing translational check)


## GTEx check

Likewise, AltNeo-B will carry out GTEx check as well, please refer to AltNeo-T documentaion for details. Also, we are able to generate a nice-looking html page to visualize expression level in GTEx of a certain splicing event through running `queryGTEx_PSI.py`. And we can assign a continuous tumor specificity score to each events by using `IW_score`.

## Overlapping with human membrane genes

For AltNeo-B, we are not interested in events that expression in cytoplasma, because there's no way for a B cell to reach to it. Instead, we focus on our search down to splicing events whose translational product would be a membrane protein. We downloaded and pre-processed all the documented human membrane protain information from Uniprot, there are two files `human_membrane_proteins.txt`, which contains all Uniprot entry name for these membrane protein. We will use Uniprot Online conversion tool, implemented as function `IDmappingACC2Ensembl` to convert Ensembl ID to Uniprot Entry. Another file is the protein sequence for all the human protein and already known isoforms `uniprot_isoform.fasta`, this will be useful when we try to align obtained peptides to known isoforms later.

After this step, basically we will have a much smaller dataframe for following steps.

## Retrieve junction sequence

This step is the same as AltNeo-T, plese check the detail there.

## Match LSVs to its potential parental transcripts

If you see this step is optional in AltNeo-T, because there we only want to control falsely discovered neoantigen candidates. Here in AltNeo-B, this is a required step, because in order to determine whether a LSV could turn into a extracellular theraputic target, we need to carefully reason what full-length peptide it will possibly generate and appraise its ability to be called a good candidate. Hence, it is the first and foremost step to recorver LSVs' parental transcripts.

The step for doing that is exactly the same as AltNeo-T, so I will refer you to AltNeo-T's documentation for details. Basically, we run three iteration to recover as many LSVs as possible, then we can use ORF finding algorithms to derive optimal ORF for each recovered transcript.

## Translational check and NMD check

However, here we add two filtering mechanisms, which is vital and highlights for AltNeo-B pipeline. Since in AltNeo-B program, we are more concerned about whether the resultant peptides derived from a LSV can finally reside on membrane, so we are not satisfied with just generating all potential ORF and ORF-translated peptides, we want to know,

1. Whether the transcript is a protein-coding trancript? If not, chances are that it might not become a valid membrane protein. (**Translational check**)

2. Whether the resultant peptides will undergo Nonsense-Mediated Decay? (**NMD**) (**NMD check**)

For the translational check, it is fairly easy, because we just need to pull out the annotation in Ensembl for each gene and see if they are protein-coding trancripts. The result will be a column like this ['', '', '#', '', '', '', '#', '', '', '', '', '', '#', ''], here each item in list means a matching endeavor to a certain transcript, emtpy string basically means it is not able to fit in a certain transcript. But for each matched item, we will label them whether we are protein-coding transcript (**#**) or not (**\***).

For NMD check, it is very involved and needs to be illustrated here. First we need to understand when NMD will happen. Basically, NMD will happen if early stop codon exists, it will prevent our body to produce too short peptides, it also server as a negative regulation mechanisms to degrade some junk peptides. To be more concret, as illustrated in **Figure 2**, Let's say we have a transcript which have 15 exons, if the stop codon resides on exon13 (**A**), which is two exon upstream of the last one (15 - 13 = 2), it will undergo NMD, because it is a early stopping. In contrast, it won't undergo NMD in situation (**B**) because it is only one exon away from last one (15 - 14 = 1). Hence, In AltNeo-B program, if stop condon-Residing Exon number is two exons away from last exon, it will labelled as (**\***), it won't affect the downstream analysis but it is a warning when you interpret the final result. However, It is noteworthy that situation is a bit more complicated than what I implemented, as shown in Situation (**C**), event though here the stop codon is on exon14, but it stops at very early part of exon14, it may still undergo NMD, but we currently don't consider that.

**Figure 2**

![Figure2](https://github.com/frankligy/AltNeo-BT/blob/main/images/B-Figure2.png)

We want to shed light on how we actually implement NMD check, it is handled by function `build_sorted_exons`, this function will take a series of subexons that constitute a certain transcript, for example, `E1.1|E1.2|E1.3|E3.1|E3.2|E4.1`, it will return a list storing all the end positions of each **"exon"**, here I need first emphasize the definition of **"exon"**, `E1.1|E1.2|E1.3` will be an exon, that's easy to understand. But in the following example, `E14.1|E14.2|E14.4|E14.5` will have two **"exons"** because there's a breaking point, a missing `E14.3`. With that understood, we now explain what `build_sorted_exons` return, as illustrated in **Figure 3**, it basically return the end point of each **"exon"**, plus the first 0. To discren this weired **"exon"** concept, we have `check_consecutive` function to take care of that which make sure the list `build_sorted_exons` function returns will consider breaking poitn situation. Having the positions, we can simply put our stop codon position onto this list, figuring out which exon this stop codon resides on, therefore to determine whether it will undergo NMD or not.

**Figure 3**

![Figure3](https://github.com/frankligy/AltNeo-BT/blob/main/images/B-Figure3.png)

After that, we check whether the junction sequence is in the range of optimal ORF we predicted out, we use variable `involvement` to represent that. In the final result, we will see a column named `involvement`, it basically mean whether a certain matching (LSV to a transcript) is valid in the sense that junction fall into the ORF range, detail for this was explained in AltNeo-T. Empty string "", which means no match, will certainly receive a False involvement.


## Alignment to known Uniprot transcripts

The unique part of AltNeo-B is that, we want to discover novel neoepitopes, we are not interested in ones that have been documented in Uniprot, which imply that they are expressed in normal tissue which prevent us from using it as a cancer immunotherapy target. Till now, we've already generated a series of possible novel isoforms and their encoded peptide sequences, but we need to verify whether they are the same as what have been documented in Uniprot dataset or not, we are interested in the latter ones.

To achieve that, we compare each derived peptides to all documented isoforms in Uniprot for a certain gene (surface protein coding gene). We adopted a seed alignment strategy, we chop each predicted peptide into 10mer and align each 10mer to each documented isoform, as shown in **Figure 4**, it will return a dictionary like this, and all predicted peptides constitute a list :

**Figure 4**

![Figure 4](https://github.com/frankligy/AltNeo-BT/blob/main/images/B-Figure4.png)

```
["Either splicing site is not involved in ORF formation or splicing event doesn't occur in this certain transcript", 
{'Q3KPI0': [False, False, False, False, 'notAligned'], 
'Q3KPI0-2': [False, False, False, False, 'notAligned'], 
'Q3KPI0-3': [False, False, False, False, 'notAligned']}, 
```


In the above example, we have two predicted peptide, one is labelled as `Either splicing site is not involved in ORF formation or splicing event doesn't occur in this certain transcript`, this means either the LSV-trancript matching is not successful, or involvement check suggested that junction won't participate in peptide formation, let alone to form neoantigen, which means we don't need to consider this event. The second peptide, however, undergo this seed alignment process to three documented isoforms, each 10mer alignment result would be a boolean variable, and the final element in the value list summarize the result for this alignment, `totallyaligned` means nothing novel, `partiallyaligned` and `notAligned` is what we are interested in. 

We will further summarize the second predicted peptide and answer the question, what kind of `identity` this predicted peptide is? Is it a `novel isofrm` that we need to carry on our pursuing or `one of already documented isoforms` that we are not interested in. As long as `totallyaligned` occur ones, it will be `one of already documented isoforms`. This result is shown in `identity` column in the final result.

## TMHMM to determine whether the aberrant region is extracellular


Now we are close but one step remains is that, the predicted peptide seems to be a good novel isoform, whether this peptide will have a extracellular domain? Because if the aberrant part is in cytoplamic side, it doesn't work either. We simply call for TMHMM algorithm to predict whether this predicted peptide contains extracellular region, this will be the final step AltNeo-B will perform for you, the result would look like this in `interest` column:

`[(True, 1), (True, 3), (True, 7), (True, 9)]`

It means check the 2nd, 4th, 8th, 10th matching and mannual double check if they are good candidates, remember Python is 0-based language. If the result is basically `False`, then no need to pursue with this splicing event.

## Novel splicing events that can not match with any existing transcripts

There are some splicing events that we are not able to match up with any exisiting transcripts, they will be sifted out when we check if all matching are empty. For those, their `involvment` column will all be `False` but we need to discern that cause them to be `False`, if they didn't pass the involvment test, that's bad, we don't want them, if they are just novel splicing type like trans-splicing event that being labelled as `False` simply because it is not able to match to existing transcripts in any way, we will further mannually check those. We use function `diffNovelFromNotInvolved_new` to achieve this goal.