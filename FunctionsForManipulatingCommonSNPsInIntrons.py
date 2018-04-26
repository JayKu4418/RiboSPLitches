from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
# from Bio import motifs
import re
from operator import itemgetter
from collections import Counter
import subprocess
import os
import sys
import cStringIO
import matplotlib.pyplot as plt
import pandas as pd

base_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}


def grabFastaSequencesForSNPWriteToFile(snpfile, window, seqdirectory, snpfastafile):
    seqtowrite = []

    chrms = range(1, 23)
    chrms.extend(['X', 'Y'])
    # print(chrms)
    chrms_act = ['chr' + str(i) for i in chrms]

    snps = pd.read_csv(snpfile,sep="\t",header=None)

    for chromosome in chrms_act:
        print(chromosome)
        snps_chrom = snps[snps[0]==chromosome]
        # print('Here')
        recsforchrom = []
        for record in SeqIO.parse(open(seqdirectory + '/' + chromosome + '.fa', 'rU'), 'fasta'):
            recsforchrom.append(record)
        print('Here2')
        for p in range(snps_chrom.shape[0]):
            snp = snps_chrom.iloc[p,:]
            if snp[4] == '-':
                seq = SeqRecord(
                    Seq(str(recsforchrom[0].seq).upper()[snp[6]-window:snp[7]+window], IUPAC.unambiguous_dna).reverse_complement(),
                    id=snp[8],
                    description=snp[8] + ', ' + chromosome + ', start = ' + str(snp[6]-window)
                                + ', end = ' + str(snp[7]+window) + ', strand = -')
            else:
                seq = SeqRecord(
                    Seq(str(recsforchrom[0].seq).upper()[snp[6] - window:snp[7] + window],
                        IUPAC.unambiguous_dna).reverse_complement(),
                    id=snp[8],
                    description=snp[8] + ', ' + chromosome + ', start = ' + str(snp[6] - window) + ', end = ' + str(
                        snp[7] + window) + ', strand = +')

            seqtowrite.append(seq)
            # print('Here3')
    SeqIO.write(seqtowrite, snpfastafile, 'fasta')


def getTrueSpliceSiteCoordsClosestToSNP(intersectingfile, end, writefile):

    # Grab all the mutations in the intersecting file with all column values
    allsnps = pd.read_csv(intersectingfile,header=None,sep="\t")

    allsnps_plus = allsnps[allsnps[4]=="+"]
    allsnps_minus = allsnps[allsnps[4]=="-"]

    if end=="5p":
        print end
        allsnps_data = {"chrom":pd.concat([allsnps_plus[0],allsnps_minus[0]]),
                           "start":pd.concat([allsnps_plus[1]-3,allsnps_minus[2]-6]),
                           "end":pd.concat([allsnps_plus[1]+6,allsnps_minus[2]+3]),
                           "snp_id":pd.concat([allsnps_plus[8],allsnps_minus[8]]),
                           "strand":pd.concat([allsnps_plus[4],allsnps_minus[4]]),
                        "score":[0]*len(pd.concat([allsnps_plus[8],allsnps_minus[8]]))}
    else:
        print end
        allsnps_data = {"chrom": pd.concat([allsnps_plus[0], allsnps_minus[0]]),
                           "start": pd.concat([allsnps_plus[2] - 20, allsnps_minus[1] - 3]),
                           "end": pd.concat([allsnps_plus[2] + 3, allsnps_minus[1] + 20]),
                           "snp_id": pd.concat([allsnps_plus[8], allsnps_minus[8]]),
                           "strand": pd.concat([allsnps_plus[4], allsnps_minus[4]]),
                        "score":[0]*len(pd.concat([allsnps_plus[8],allsnps_minus[8]]))}

    allsnps_df = pd.DataFrame(allsnps_data,columns=["chrom","start","end","snp_id","score","strand"])

    allsnps_df.to_csv(writefile,header=False,index=False,sep="\t")

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


# This function will return mutations whose mutated sequence has a maxentscore that is comparable to that of the true splice sites
# Remember for calculating max ent score of 5' splice site, you need last three bases of exon and first six bases of intron
# For 3' splice site, you need last 20 bases of intron and first three bases of exon
def returnSNPsCreatingSpliceSitesWithComparableScoresAsTrueSpliceSite(snpseqfile, snpsfile, truesplicesitefile,
                                                                        splicesite, min_sites,
                                                                           comparablescore):
    # Grab all the records from mutseqfile store in variable allrecs
    allrecs = []
    for record in SeqIO.parse(open(snpseqfile, 'rU'), 'fasta'):
        allrecs.append(record)

    # Grab all the mutations from the mutfile
    allsnps = pd.read_csv(snpsfile,header=None,sep="\t")

    # Grab all the true splice sites store in var list truesplicesites
    truesplicesites = []
    for record in SeqIO.parse(open(truesplicesitefile, 'rU'), 'fasta'):
        truesplicesites.append(record)

    recs_to_consider = []
    # Go through every record in allrecs
    for rec in allrecs:
        # This is the list that will store the various arrangments of sequences for that rec
        seqs_for_rec = []
        # Grab the id and sequence for the rec
        rec_snpname = rec.id
        print rec_snpname
        # Since sequence around mutations were stored with a 50 bp window, we only need the 20 base pair window for
        # 5' splice sites
        if splicesite == '5p':
            rec_seq = rec.upper().seq[15:36]
        else:
            rec_seq = rec.upper().seq

        # Find the mutation that matches the rec and if no match is found, break out of loop
        snps_for_rec = allsnps[allsnps[8]==rec_snpname].values.tolist()
        print snps_for_rec
        if len(snps_for_rec) == 0:
            print("No snp found for record")
            print rec_snpname
            break
        else:
            for z in range(len(snps_for_rec)):
                snp_for_rec = snps_for_rec[z]
                # If strand is minus, we need to get complement of base and it's mutation, if not just keep base.
                # Use dictionary to grab complement
                if snp_for_rec[4] == "-":
                    mutated_base = base_dict[snp_for_rec[9]]
                    base_mutated_to = base_dict[snp_for_rec[10]]
                    print mutated_base
                    print base_mutated_to
                else:
                    mutated_base = snp_for_rec[9]
                    base_mutated_to = snp_for_rec[10]
                # If the base of concern on the sequence of the record is not the same as the mutated base, break out of loop
                if (splicesite == "5p" and rec_seq[10] != mutated_base) or (
                        splicesite == "3p" and rec_seq[25] != mutated_base):
                    print("WT mutation site does not match record")
                    print rec_snpname
                    break
                else:
                    # Have four different 9 base sequences, where base of concern is either last 2 bases of exon or first two bases
                    # of intron, and the corresponding sequence with the mutation.
                    if splicesite == "5p":
                        if min_sites:
                            for i in range(0, 4):
                                # When i = 0, the mutation site is at last base of the first two bases of the intron,
                                # as i increases, the site will then get closer to the first base of the last two bases
                                # of the exon
                                seq_i = rec_seq[6 + i:15 + i]
                                seq_mutated_i = seq_i[0:4 - i] + base_mutated_to + seq_i[5 - i:]
                                seqs_for_rec.append(SeqRecord(seq_i, id=rec_snpname + "_" + str(i + 1) + "_WT"))
                                seqs_for_rec.append(SeqRecord(seq_mutated_i, id=rec_snpname + "_" + str(i + 1) + "_MUT"))
                        else:
                            for i in range(0, 9):
                                # When i = 0, the mutation site is at last base of the first six bases of the intron,
                                # as i increases, the site will then get closer to the first base of the last two bases
                                # of the exon
                                seq_i = rec_seq[2 + i:11 + i]
                                if i == 0:
                                    seq_mutated_i = seq_i[0:8] + base_mutated_to
                                elif i == 8:
                                    seq_mutated_i = base_mutated_to + seq_i[1:]
                                else:
                                    seq_mutated_i = seq_i[0:8 - i] + base_mutated_to + seq_i[9 - i:]
                                seqs_for_rec.append(SeqRecord(seq_i, id=rec_snpname + "_" + str(i + 1) + "_WT"))
                                seqs_for_rec.append(SeqRecord(seq_mutated_i, id=rec_snpname + "_" + str(i + 1) + "_MUT"))
                    else:
                        if min_sites:
                            for i in range(0, 4):
                                # When i = 0, the mutation site is at last base of the first two bases of the intron,
                                # as i increases, the site will then get closer to the first base of the last two bases
                                # of the exon
                                seq_i = rec_seq[4 + i:27 + i]
                                seq_mutated_i = seq_i[0:21 - i] + base_mutated_to + seq_i[22 - i:]
                                seqs_for_rec.append(SeqRecord(seq_i, id=rec_snpname + "_" + str(i + 1) + "_WT"))
                                seqs_for_rec.append(SeqRecord(seq_mutated_i, id=rec_snpname + "_" + str(i + 1) + "_MUT"))
                        else:
                            for i in range(0, 23):
                                # When i = 0, the mutation site is at last base of the first three bases of the exon,
                                # as i increases, the site will then get closer to the first base of the last 20 bases
                                # of the intron
                                seq_i = rec_seq[3 + i:26 + i]
                                if i == 0:
                                    seq_mutated_i = seq_i[0:22] + base_mutated_to
                                elif i == 22:
                                    seq_mutated_i = base_mutated_to + seq_i[1:]
                                else:
                                    seq_mutated_i = seq_i[0:22 - i] + base_mutated_to + seq_i[23 - i:]
                                seqs_for_rec.append(SeqRecord(seq_i, id=rec_snpname + "_" + str(i + 1) + "_WT"))
                                seqs_for_rec.append(SeqRecord(seq_mutated_i, id=rec_snpname + "_" + str(i + 1) + "_MUT"))

                    seqs_for_rec.extend([i for i in truesplicesites if i.id.split('_')[0] == rec_snpname])
                    # Write these 8 sequences into a fasta file
                    #SeqIO.write(seqs_for_rec, '../../MaxEntScore/RecSeqsToAnalyze.fa', 'fasta')
                    SeqIO.write(seqs_for_rec, 'RecSeqsToAnalyze.fa', 'fasta')
                    # Now run the perl script for MaxEntScore
                    if splicesite == "5p":
                        pipe_out = subprocess.check_output(['score5.pl', 'RecSeqsToAnalyze.fa'])
                    else:

                        pipe_out = subprocess.check_output(['score3.pl', 'RecSeqsToAnalyze.fa'])

                    #os.chdir('../Projects/Generating_Statistics_About_Introns')

                    maxentscores_for_rec_mut_pos = []

                    # Grab the max ent scores for the sequences and only keep the mutation if max ent score of the sequences
                    # with mutation are greater than zero
                    maxentscores_for_rec = [i.split('\t') for i in pipe_out.strip().split('\n')]
                    if min_sites:
                        maxentscores_for_rec_mut = [maxentscores_for_rec[i] for i in range(1, 8, 2)]
                        for j in maxentscores_for_rec[8:]:
                            maxentscores_for_rec_mut_pos.extend(
                                [i[1] for i in maxentscores_for_rec_mut if float(i[1]) >= float(j[1]) - comparablescore])
                    else:
                        if splicesite == "5p":
                            maxentscores_for_rec_mut = [maxentscores_for_rec[i] for i in range(1, 18, 2)]
                            for j in maxentscores_for_rec[18:]:
                                maxentscores_for_rec_mut_pos.extend([i[1] for i in maxentscores_for_rec_mut if
                                                                     float(i[1]) > 0 and
                                                                     float(i[1]) >= float(j[1]) - comparablescore])
                        else:
                            maxentscores_for_rec_mut = [maxentscores_for_rec[i] for i in range(1, 46, 2)]
                            for j in maxentscores_for_rec[46:]:
                                maxentscores_for_rec_mut_pos.extend([i[1] for i in maxentscores_for_rec_mut if
                                                                     float(i[1]) > 0 and
                                                                     float(i[1]) >= float(j[1]) - comparablescore])

                    if len(maxentscores_for_rec_mut_pos) > 0:
                        recs_to_consider.append(rec_snpname)
    #if os.getcwd() != '/Users/jayashreekumar/Documents/LaederachLab/Projects/Generating_Statistics_About_Introns':
        #os.chdir('../Projects/Generating_Statistics_About_Introns')
    return recs_to_consider
