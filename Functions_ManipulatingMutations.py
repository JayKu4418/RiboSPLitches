from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
# from Bio import motifs
import regex as re2
from operator import itemgetter
from collections import Counter
import subprocess
import os
import sys
import cStringIO
import matplotlib.pyplot as plt
import pandas as pd

base_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

def getTrueSpliceSiteCoordsClosestToMut(intersectingfile, end, writefile):

    # Grab all the mutations in the intersecting file with all column values
    allmuts = pd.read_csv(intersectingfile,header=None,sep="\t")

    allmust_plus = allmuts[allmuts[4]=="+"]
    allmust_minus = allmuts[allmuts[4]=="-"]

    if end=="5p":
        print end
        allmuts_data = {"chrom":pd.concat([allmust_plus[0],allmust_minus[0]]),
                           "start":pd.concat([allmust_plus[1]-3,allmust_minus[2]-6]),
                           "end":pd.concat([allmust_plus[1]+6,allmust_minus[2]+3]),
                           "snp_id":pd.concat([allmust_plus[9],allmust_minus[9]]),
                           "strand":pd.concat([allmust_plus[4],allmust_minus[4]]),
                        "score":[0]*len(pd.concat([allmust_plus[9],allmust_minus[9]]))}
    else:
        print end
        allmuts_data = {"chrom": pd.concat([allmust_plus[0], allmust_minus[0]]),
                           "start": pd.concat([allmust_plus[2] - 20, allmust_minus[1] - 3]),
                           "end": pd.concat([allmust_plus[2] + 3, allmust_minus[1] + 20]),
                           "snp_id": pd.concat([allmust_plus[9], allmust_minus[9]]),
                           "strand": pd.concat([allmust_plus[4], allmust_minus[4]]),
                        "score":[0]*len(pd.concat([allmust_plus[9],allmust_minus[9]]))}

    allmuts_df = pd.DataFrame(allmuts_data,columns=["chrom","start","end","snp_id","score","strand"])

    allmuts_df.to_csv(writefile,header=False,index=False,sep="\t")

def getMutationsThatCreateOrDestroyISREMotifs(ISREfile, mutation_seqfile, mutationsOfInterest, mutsite, create):
    # Open the file containing ISRE motifs
    ISRE_motifs = pd.read_csv(ISREfile,header=None)
    print ISRE_motifs.head()
    # Open file with mutations and sequences surrounding mutations
    mutations_withSeq = pd.read_csv(mutation_seqfile,sep="\t",header=None)
    print mutations_withSeq.shape
    # Open file that contains IDs of mutations of interest
    muts_Of_interest = pd.read_csv(mutationsOfInterest,sep="\t",header=None)
    print muts_Of_interest.shape
    
    # Subset muations_withSeq data and only have ones that contains IDs of mutations of interest
    muts_Of_interest_withSeq = mutations_withSeq[mutations_withSeq[9].isin(muts_Of_interest[9])]
    print muts_Of_interest_withSeq.shape
    
    if create == True:
        index_for_seq = 14
    else:
        index_for_seq = 13
    
    # Create a dictionary that stores the mutated sequence as values with mutation id as the key
    muts_dict = pd.Series(muts_Of_interest_withSeq[index_for_seq].values,index=muts_Of_interest_withSeq[9].values).to_dict()
    
    # These are the lists that are going to contain the mutation ids for mutations that result in new ISRE motifs
    mutations_Have_Motifs = []

    # Go through every record in mutation file
    for mut in muts_dict.keys():
        mut_seq = muts_dict[mut]
        mutatedISREmotifs_for_mut = []
        for i in ISRE_motifs[0].values:
            x = [[m.start(), m.end()] for m in re2.finditer(i, mut_seq,overlapped=True) if mutsite in range(m.start(), m.end())]
            if len(x) != 0:
                mutatedISREmotifs_for_mut.append(i)
        if len(mutatedISREmotifs_for_mut) > 0:
            mutations_Have_Motifs.append(mut)

    return mutations_Have_Motifs