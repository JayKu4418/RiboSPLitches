# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 12:48:56 2016

@author: jayashreekumar
"""

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


def grabFastaSequencesForMutWriteToFile(mutfile, window, seqdirectory, mutfastafile):
    seqtowrite = []

    chrms = range(1, 23)
    chrms.extend(['X', 'Y'])
    # print(chrms)
    chrms_act = ['chr' + str(i) for i in chrms]

    for chromosome in chrms_act:
        print(chromosome)
        with open(mutfile) as f:
            mutstograb = [[i.strip().split('\t')[0], int(i.strip().split('\t')[2]) - window,
                           int(i.strip().split('\t')[2]) + window, i.strip().split('\t')[3],
                           i.strip().split('\t')[8].split(';')[3].split('=')[1]] for i in f if
                          i.strip().split('\t')[0] == chromosome]
        # print('Here')
        recsforchrom = []
        for record in SeqIO.parse(open(seqdirectory + '/' + chromosome + '.fa', 'rU'), 'fasta'):
            recsforchrom.append(record)
        print('Here2')
        for p in range(len(mutstograb)):
            mut = mutstograb[p]
            if mut[4] == '-':
                seq = SeqRecord(
                    Seq(str(recsforchrom[0].seq).upper()[mut[1] - 1:mut[2]], IUPAC.unambiguous_dna).reverse_complement(),
                    id=mut[3],
                    description=mut[3] + ', ' + chromosome + ', start = ' + str(mut[1]) + ', end = ' + str(mut[2]))
            else:
                seq = SeqRecord(Seq(str(recsforchrom[0].seq).upper()[mut[1] - 1:mut[2]], IUPAC.unambiguous_dna), id=mut[3],
                                description=mut[3] + ', ' + chromosome + ', start = ' + str(mut[1]) + ', end = ' + str(
                                    mut[2]))

            seqtowrite.append(seq)
            # print('Here3')
    SeqIO.write(seqtowrite, mutfastafile, 'fasta')


def grabFastaSequencesForBedFileWriteToFile(bedfile, seqdirectory, bedfastafile):
    seqtowrite = []

    chrms = range(1, 23)
    chrms.extend(['X', 'Y'])
    chrms_act = ['chr' + str(i) for i in chrms]

    for chromosome in chrms_act:
        print(chromosome)
        with open(bedfile) as f:
            bedlines = [[i.strip().split('\t')[0], int(i.strip().split('\t')[1]), int(i.strip().split('\t')[2]),
                         i.strip().split('\t')[3]] for i in f if i.strip().split('\t')[0] == chromosome]
        recsforchrom = []
        for record in SeqIO.parse(open(seqdirectory + '/' + chromosome + '.fa', 'rU'), 'fasta'):
            recsforchrom.append(record)
        print('Here2')
        for p in range(len(bedlines)):
            x = bedlines[p]
            x_seq = Seq(str(recsforchrom[0].seq).upper()[x[1] - 1:x[2]], IUPAC.unambiguous_dna)
            if x[3] == '+':
                seqToWrite = x_seq
            else:
                seqToWrite = x_seq.reverse_complement()

            seq = SeqRecord(seqToWrite, id=chromosome + '_' + str(p + 1),
                            description=chromosome + ', start= ' + str(x[1]) + ', end= ' + str(
                                x[2]) + ', strand= ' + str(x[3]))
            seqtowrite.append(seq)

            # print('Here3')
    SeqIO.write(seqtowrite, bedfastafile, 'fasta')


def findSeqInFasta(seqtofind, window, fastafile, writefile):
    allseq = []

    for record in SeqIO.parse(open(fastafile, 'rU'), 'fasta'):
        rec_seq = str(record.upper().seq)
        if seqtofind in rec_seq:
            seq_index_rec = [(m.start(), m.end()) for m in re.finditer(seqtofind, rec_seq)]
            seq_for_rec = [rec_seq[(i[0] - window):(i[1] + window)] for i in seq_index_rec if
                           (i[0] - window) >= 0 and (i[0] + window) < len(rec_seq)]
            allseq.extend(seq_for_rec)

    allseq_seq = [SeqRecord(Seq(allseq[i], IUPAC.unambiguous_dna), id=seqtofind + '_' + str(i)) for i in
                  range(len(allseq))]

    SeqIO.write(allseq_seq, writefile, 'fasta')


def grabUnique_NonOverlapping_Coords(bedfile, writefile):
    with open(bedfile) as f:
        all_lines = [line.strip().split('\t') for line in f]

    chrms = range(1, 23)
    chrms.extend(['X', 'Y'])
    chrms_act = ['chr' + str(i) for i in chrms]

    with open(writefile, 'w') as fw:
        for chromosome in chrms_act:
            print(chromosome)
            # Grab lines for that chromosome from bedfile
            lines_chr = [i for i in all_lines if i[0] == chromosome]

            # This vector will contain  all the introns that can be contained within another intron
            trial_1 = []

            # Go through every intron in bed file
            for i in lines_chr:
                # Collect all lines that are encompassed within current line i
                x = [j for j in lines_chr if
                     j != i and j[3] == i[3] and int(j[1]) >= int(i[1]) and int(j[2]) <= int(i[2])]
                # add those lines which are encompassed by another line i
                if len(x) > 0:
                    trial_1.extend(x)

            # Grab only lines that are not included in trial_1, so these will be the introns that encompass others.
            # Sort these and divide between plus and minus strand
            lines_to_include_chr = [[i[0], int(i[1]), int(i[2]), i[3]] for i in lines_chr if i not in trial_1]
            lines_to_include_chr_sorted = sorted(lines_to_include_chr, key=itemgetter(1))
            lines_to_include_chr_sorted_plus = [i for i in lines_to_include_chr_sorted if i[3] == '+']
            lines_to_include_chr_sorted_minus = [i for i in lines_to_include_chr_sorted if i[3] == '-']

            # we need to now join overlapping introns together as one
            # Create an empty vector to store the final set of lines
            final_set_lines_chr = []

            # We first look to see if the end coordinate of the current line is within the start and end coordinate of the next line
            # If it is, we keep looking at the following lines to see if the previous lines end coordinate is within the start and end of the next line
            # When this is no longer the case, we create one long intron where the start coordinate is from the first line and the end coordinate is from the final line
            # If end of coordinate of current line is not within the start and end coordinate of next line, then we just add current line to our final list

            # This is for the plus strand
            j = 0
            while j < len(lines_to_include_chr_sorted_plus):
                if j != len(lines_to_include_chr_sorted_plus) - 1 and lines_to_include_chr_sorted_plus[j][2] >= \
                        lines_to_include_chr_sorted_plus[j + 1][1] and lines_to_include_chr_sorted_plus[j][2] <= \
                        lines_to_include_chr_sorted_plus[j + 1][2]:

                    i = j
                    while i != len(lines_to_include_chr_sorted_plus) - 1 and lines_to_include_chr_sorted_plus[i][2] >= \
                            lines_to_include_chr_sorted_plus[i + 1][1] and lines_to_include_chr_sorted_plus[i][2] <= \
                            lines_to_include_chr_sorted_plus[i + 1][2]:
                        i = i + 1

                    newline = [lines_to_include_chr_sorted_plus[j][0], lines_to_include_chr_sorted_plus[j][1],
                               lines_to_include_chr_sorted_plus[i][2], '+']

                    final_set_lines_chr.append(newline)
                    j = i + 1
                else:
                    final_set_lines_chr.append(lines_to_include_chr_sorted_plus[j])
                    j = j + 1

            # This is for the plus strand
            j = 0
            while j < len(lines_to_include_chr_sorted_minus):
                if j != len(lines_to_include_chr_sorted_minus) - 1 and lines_to_include_chr_sorted_minus[j][2] >= \
                        lines_to_include_chr_sorted_minus[j + 1][1] and lines_to_include_chr_sorted_minus[j][2] <= \
                        lines_to_include_chr_sorted_minus[j + 1][2]:

                    i = j
                    while i != len(lines_to_include_chr_sorted_minus) - 1 and lines_to_include_chr_sorted_minus[i][2] >= \
                            lines_to_include_chr_sorted_minus[i + 1][1] and lines_to_include_chr_sorted_minus[i][2] <= \
                            lines_to_include_chr_sorted_minus[i + 1][2]:
                        i = i + 1

                    newline = [lines_to_include_chr_sorted_minus[j][0], lines_to_include_chr_sorted_minus[j][1],
                               lines_to_include_chr_sorted_minus[i][2], '+']

                    final_set_lines_chr.append(newline)
                    j = i + 1
                else:
                    final_set_lines_chr.append(lines_to_include_chr_sorted_minus[j])
                    j = j + 1

            # Sort the final set of lines and write to file
            final_set_lines_chr_sorted = sorted(final_set_lines_chr, key=itemgetter(1))
            for i in final_set_lines_chr_sorted:
                fw.write('\t'.join([i[0], str(i[1]), str(i[2]), '-', i[3]]))
                fw.write('\n')


# This function will write mutations that are flanking ISRE motifs within a certain distance to a file
# and return those don't flank motifs in a list
def returnMutationsFlankingMotif(ISREfile, mutation_seqfile, mutation_file, flankdist,
                                 mutation_flankingsmutfile):
    # Open the file containing ISRE motifs
    with open(ISREfile) as f:
        ISRE_motifs = [line.strip().split('\t')[0:2] for line in f]

    # These lists will contain mutations that flank motifs and those that don't
    mutsites_flanking_motifs = []
    mutsite_noflanking_motifs = []

    # This will allow us to figure out the right index for the mutation site
    # sites have 51 bases with index 25 being the mut site
    mutsite = 25

    # Go through every record in file
    for record in SeqIO.parse(open(mutation_seqfile, 'rU'), 'fasta'):
        # print record.seq
        rec_up = record.upper()
        flanking_motifs_in_record = []
        for i in ISRE_motifs:
            # print i
            # print record.seq
            x = [[m.start(), m.end() - 1] for m in re.finditer(i[0], str(rec_up.seq))]
            # print x
            # print len(x)
            # If motif is found in the sequence, calculate how far the mutation site is from the motif, both from
            # the start of the motif and end of the motif and then grab the minimum of that
            if len(x) != 0:
                dist_flanking_for_i = []
                for j in x:
                    mut_dist_from_motif = min([abs(j[0] - mutsite), abs(j[1] - mutsite)])
                    if mut_dist_from_motif <= flankdist:
                        if abs(j[0] - mutsite) < abs(j[1] - mutsite):
                            dist_flanking_for_i.append(-mut_dist_from_motif)
                        else:
                            dist_flanking_for_i.append(5+mut_dist_from_motif)
                if len(dist_flanking_for_i)!=0:
                    flanking_motifs_in_record.append([i,dist_flanking_for_i])
        if len(flanking_motifs_in_record) != 0:
            mutsites_flanking_motifs.append([record.id, rec_up.seq, flanking_motifs_in_record])
        # If mutation site does not occur within flanking of motifs
        else:
            mutsite_noflanking_motifs.append([record.id, rec_up.seq])

    with open(mutation_file) as f:
        muts = [line.strip().split('\t') for line in f]

    allflankingmuts_ids = [i[0] for i in mutsites_flanking_motifs]
    muts_to_write = [i for i in muts if i[3] in allflankingmuts_ids]

    if mutation_flankingsmutfile != "":
        with open(mutation_flankingsmutfile, 'w') as fw:
            for j in muts_to_write:
                fw.write('\t'.join(j))
                fw.write('\n')
    return [mutsites_flanking_motifs, mutsite_noflanking_motifs]


# Separate mutation sites found within unmutated motifs and mutations that create motifs, return those in a list
# The mutation sites that are not found in any motifs mutated or unmutated, write those into a file
def separateMutationSites(ISREfile, mutation_seqfile, mutationfile, strand,
                          mutation_possible_flanking_seqfile, mutation_possible_flanking_file):
    # Open the file containing ISRE motifs
    with open(ISREfile) as f:
        ISRE_motifs = [line.strip().split('\t')[0:2] for line in f]
    with open(mutationfile) as f:
        mutations = [line.strip().split('\t') for line in f]

    # Create a dictionary that stores the WT and MUT base as values with mutation id as the key
    muts_dict = {}
    for i in mutations:
        muts_dict[i[3]] = [i[5], i[6]]

    # These are the lists that are going to contain the mutation ids for unmutated, mutated and noncontaining motifs
    mutsites_contain_unmutated_motifs = []
    mutsites_contain_mutated_motifs = []
    mutsites_notcontaining_motifs = []

    # This will allow us to figure out the right index for the mutation site
    # because donor sites have 21 bases with index 10 being the mut site
    # while acceptor sites have 51 bases with index 25 being the mut site
    mutsite = 25

    # Go through every record in mutation file
    for record in SeqIO.parse(open(mutation_seqfile, 'rU'), 'fasta'):
        # print record.seq
        rec_up = record.upper()
        # Looking if mutation site occurs within ISRE motif
        motifs_in_record_unmutated = []
        # Look at every motif and see if the mutated site is found within an existing motif
        for i in ISRE_motifs:
            x = [[m.start(), m.end()] for m in re.finditer(i[0], str(rec_up.seq)) if
                 mutsite in range(m.start(), m.end())]
            if len(x) != 0:
                dist_from_motif = [25-p[0] for p in x]
                motifs_in_record_unmutated.append([i,dist_from_motif])
        if len(motifs_in_record_unmutated) != 0:
            mutsites_contain_unmutated_motifs.append([record.id, rec_up.seq, motifs_in_record_unmutated])
        # If mutation site does not occur within ISRE motif
        else:
            # Look if you use the mutated base, does the mutation site occur within ISRE motif
            # Change original base to mutated base
            mut_change = muts_dict[record.id]
            if strand == "minus":
                mutant = [base_dict[mut_change[0]], base_dict[mut_change[1]]]
            else:
                mutant = mut_change
            seq = str(rec_up.seq)
            if seq[mutsite] == mutant[0]:
                seqtomutate = seq[0:mutsite] + mutant[1] + seq[mutsite + 1:]
                motifs_in_record_mutated = []
                for i in ISRE_motifs:
                    x = [[m.start(), m.end()] for m in re.finditer(i[0], seqtomutate) if mutsite in range(m.start(), m.end())]
                    if len(x) != 0:
                        dist_from_motif = [25 - p[0] for p in x]
                        motifs_in_record_mutated.append([i,dist_from_motif])
                if len(motifs_in_record_mutated) != 0:
                    mutsites_contain_mutated_motifs.append([record.id, rec_up.seq, motifs_in_record_mutated])
                # If mutated base does not occur within motif
                else:
                    # flanking_motifs_in_record = []
                    # for i in ISRE_motifs:
                    #     #print i
                    #     #print record.seq
                    #     x = [m for m in re.finditer(i[0],str(rec_up.seq)) if 10 in range(m.start()-3,m.end()+3)]
                    #     if len(x)!=0:
                    #         flanking_motifs_in_record.append(i)
                    # if len(flanking_motifs_in_record)!=0:
                    #     mutsites_flanking_motifs.append([record.id,rec_up.seq,flanking_motifs_in_record])
                    # # If mutation site does not occur within flanking of motifs
                    # else:
                    #     mutsites_no_motifs.append([record.id,rec_up.seq])
                    mutsites_notcontaining_motifs.append(record)
            else:
                print "WT base in sequence did not match the WT base in the mutation record"
                print mutant
                print(record.id)

    SeqIO.write(mutsites_notcontaining_motifs, mutation_possible_flanking_seqfile, 'fasta')

    records_to_write = [rec.id for rec in mutsites_notcontaining_motifs]

    mutations_to_write = [i for i in mutations if i[3] in records_to_write]

    with open(mutation_possible_flanking_file, 'w') as fw:
        for j in mutations_to_write:
            fw.write(('\t').join(j))
            fw.write('\n')

    return [mutsites_contain_unmutated_motifs, mutsites_contain_mutated_motifs,[k.id for k in mutsites_notcontaining_motifs]]

def getMutationsThatCreateISREMotifs(ISREfile, mutation_seqfile, mutationfile, mutation_possible_flanking_file):
    # Open the file containing ISRE motifs
    with open(ISREfile) as f:
        ISRE_motifs = [line.strip().split('\t')[0:2] for line in f]
    with open(mutationfile) as f:
        mutations = [line.strip().split('\t') for line in f]

    # Create a dictionary that stores the WT and MUT base as values with mutation id as the key
    muts_dict = {}
    for i in mutations:
        muts_dict[i[3]] = [i[5], i[6],i[8].split(';')[3].split('=')[1]]

    # These are the lists that are going to contain the mutation ids for mutated and noncontaining mutated motifs
    mutsites_contain_mutated_motifs = []
    mutsites_notcontaining_mutated_motifs = []
    # This will allow us to figure out the right index for the mutation site
    # because donor sites have 21 bases with index 10 being the mut site
    # while acceptor sites have 51 bases with index 25 being the mut site
    mutsite = 25

    # Go through every record in mutation file
    for record in SeqIO.parse(open(mutation_seqfile, 'rU'), 'fasta'):
        if record.id in muts_dict.keys():
            # print record.seq
            rec_up = record.upper()
             # Look if you use the mutated base, does the mutation site occur within ISRE motif
            # Change original base to mutated base
            mut_change = muts_dict[record.id]
            if mut_change[2] == "-":
                mutant = [base_dict[mut_change[0]], base_dict[mut_change[1]]]
            else:
                mutant = mut_change[0:2]
            seq = str(rec_up.seq)
            if seq[mutsite] == mutant[0]:
                seqtomutate = seq[0:mutsite] + mutant[1] + seq[mutsite + 1:]
                motifs_in_record_mutated = []
                for i in ISRE_motifs:
                    x = [[m.start(), m.end()] for m in re.finditer(i[0], seqtomutate) if mutsite in range(m.start(), m.end())]
                    if len(x) != 0:
                        dist_from_motif = [25 - p[0] for p in x]
                        motifs_in_record_mutated.append([i,dist_from_motif])
                if len(motifs_in_record_mutated) != 0:
                    mutsites_contain_mutated_motifs.append([record.id, rec_up.seq, motifs_in_record_mutated])
                # If mutated base does not occur within motif
                else:
                    mutsites_notcontaining_mutated_motifs.append(record)
            else:
                print "WT base in sequence did not match the WT base in the mutation record"
                print mutant
                print(record.id)

    #SeqIO.write(mutsites_notcontaining_mutated_motifs, mutation_possible_flanking_seqfile, 'fasta')

    records_to_write = [rec.id for rec in mutsites_notcontaining_mutated_motifs]

    mutations_to_write = [i for i in mutations if i[3] in records_to_write]

    with open(mutation_possible_flanking_file, 'w') as fw:
        for j in mutations_to_write:
            fw.write(('\t').join(j))
            fw.write('\n')

    return [ mutsites_contain_mutated_motifs,[k.id for k in mutsites_notcontaining_mutated_motifs]]

# This is the function
def findMotifsWithinMutationAndFlanks(ISREfile, mutation_flanking_seqfile):
    # Open the file containing ISRE motifs
    with open(ISREfile) as f:
        ISRE_motifs = [line.strip().split('\t')[0:2] for line in f]

    # Find mutations that contain ISRE motifs
    muts_contains_motifs = []
    for record in SeqIO.parse(open(mutation_flanking_seqfile, 'rU'), 'fasta'):
        motifs_in_mut = [i for i in ISRE_motifs if i[1] in record.seq]
        if len(motifs_in_mut) > 0:
            muts_contains_motifs.append([record.id, record.seq, motifs_in_mut])

    print(len(muts_contains_motifs))
    # We want to separate splice mutations that happen within the motif and those flanking motifs
    recsForMut_within_motif = []
    recsForMut_not_within_motif = []
    for i in muts_contains_motifs:
        motif_for_rec = []
        motif_not_for_rec = []
        for j in i[2]:
            x = [(m.start(), m.end()) for m in re.finditer(j[1], str(i[1])) if 10 in range(m.start(), m.end())]
            if len(x) > 0:
                motif_for_rec.append(j)
            else:
                y = [(m.start(), m.end()) for m in re.finditer(j[1], str(i[1])) if
                     10 in range(m.start() - 3, m.end() + 3)]
                if len(y) > 0:
                    motif_not_for_rec.append(j)
        if len(motif_for_rec) > 0:
            recsForMut_within_motif.append([i[0], i[1], motif_for_rec])
        elif len(motif_not_for_rec) > 0:
            recsForMut_not_within_motif.append([i[0], i[1], motif_not_for_rec])

    return ([recsForMut_within_motif, recsForMut_not_within_motif])


def findMutsFlankingMotifsWithinChildrenMotifs(mutsOfInterest, ISREfile):
    muts_flankingMotifs_InChildrenMotifs = []

    with open(ISREfile) as f:
        ISREs = [line.strip().split('\t') for line in f]
    ISRE_dict = {}
    for i in ISREs:
        ISRE_dict[i[0]] = i[1:]

    for i in mutsOfInterest:
        seq_i = str(i[1])
        mut_motif_under_parent = []
        for j in i[2]:
            ISRE_j = ISRE_dict[j[0]]
            motif_child = ISRE_j[1].split(',')
            motif_child_in_mut = []
            for k in motif_child:
                if k != ' ':
                    y = [(m.start(), m.end()) for m in re.finditer(k, seq_i) if 10 in range(m.start(), m.end())]
                    if len(y) != 0:
                        motif_child_in_mut.append(k)
            if len(motif_child_in_mut) != 0:
                mut_motif_under_parent.append(motif_child_in_mut)
        if len(mut_motif_under_parent) != 0:
            muts_flankingMotifs_InChildrenMotifs.append([i[0], i[1], mut_motif_under_parent])

    return muts_flankingMotifs_InChildrenMotifs


def appendMutationsToBigList(splicemutations, bigsplicemutationfile, strand, ud):
    for i in splicemutations:
        x = [i[0], str(i[1])]
        for j in i[2]:
            x.append(j[0])
            x.append(j[1])
        x.append(strand)
        x.append(ud)
        bigsplicemutationfile.append(x)


def captureMaxEntScoreForSequenceMutations(mutseqfile, mutfile, truesplicesitefile, mut, minus, splicesite, min_sites):
    # Grab all the records from mutseqfile store in variable allrecs
    allrecs = []
    for record in SeqIO.parse(open(mutseqfile, 'rU'), 'fasta'):
        allrecs.append(record)

    # Grab all the mutations from the mutfile
    with open(mutfile) as f:
        allmuts = [line.strip().split('\t') for line in f]

    # Grab all the true splice sites store in var list truesplicesites
    truesplicesites = []
    for record in SeqIO.parse(open(truesplicesitefile, 'rU'), 'fasta'):
        truesplicesites.append(record)

    rec = [i for i in allrecs if i.id == mut][0]

    # This is the list that will store the various arrangments of sequences for that rec
    seqs_for_rec = []
    # Grab the id and sequence for the rec
    rec_mutname = rec.id
    rec_seq = rec.upper().seq

    if splicesite == '5p':
        rec_seq = rec.upper().seq[15:36]
    else:
        rec_seq = rec.upper().seq

    # Find the mutation that matches the rec
    mut_for_rec = [i for i in allmuts if i[3] == rec_mutname]
    if len(mut_for_rec) != 1:
        print("No mutation found for record")
        print rec_mutname
    else:
        # If strand is minus, we need to get complement of base and it's mutation, if not just keep base.
        # Use dictionary to grab complement
        if minus:
            mutated_base = base_dict[mut_for_rec[0][5]]
            base_mutated_to = base_dict[mut_for_rec[0][6]]
        else:
            mutated_base = mut_for_rec[0][5]
            base_mutated_to = mut_for_rec[0][6]
        # If the base of concern on the sequence of the record is not the same as the mutated base, break out of loop
        if (splicesite == "5p" and rec_seq[10] != mutated_base) or (splicesite == "3p" and rec_seq[25] != mutated_base):
            print("WT mutation site does not match record")
            print rec_mutname
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
                        seqs_for_rec.append(SeqRecord(seq_i, id=rec_mutname + "_" + str(i + 1) + "_WT"))
                        seqs_for_rec.append(SeqRecord(seq_mutated_i, id=rec_mutname + "_" + str(i + 1) + "_MUT"))
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
                        seqs_for_rec.append(SeqRecord(seq_i, id=rec_mutname + "_" + str(i + 1) + "_WT"))
                        seqs_for_rec.append(SeqRecord(seq_mutated_i, id=rec_mutname + "_" + str(i + 1) + "_MUT"))
            else:
                if min_sites:
                    for i in range(0, 4):
                        # When i = 0, the mutation site is at last base of the first two bases of the intron,
                        # as i increases, the site will then get closer to the first base of the last two bases
                        # of the exon
                        seq_i = rec_seq[4 + i:27 + i]
                        seq_mutated_i = seq_i[0:21 - i] + base_mutated_to + seq_i[22 - i:]
                        seqs_for_rec.append(SeqRecord(seq_i, id=rec_mutname + "_" + str(i + 1) + "_WT"))
                        seqs_for_rec.append(SeqRecord(seq_mutated_i, id=rec_mutname + "_" + str(i + 1) + "_MUT"))
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
                        seqs_for_rec.append(SeqRecord(seq_i, id=rec_mutname + "_" + str(i + 1) + "_WT"))
                        seqs_for_rec.append(SeqRecord(seq_mutated_i, id=rec_mutname + "_" + str(i + 1) + "_MUT"))

            seqs_for_rec.extend([i for i in truesplicesites if i.id.split('_')[0] == rec_mutname])
            # Write these 8 sequences into a fasta file
            SeqIO.write(seqs_for_rec, '../MaxEntScore/RecSeqsToAnalyze.fa', 'fasta')

            # Now run the perl script for MaxEntScore
            os.chdir('../MaxEntScore')

            if splicesite == '5p':
                pipe_out = subprocess.check_output(['perl', 'score5.pl', 'RecSeqsToAnalyze.fa'])
            else:
                pipe_out = subprocess.check_output(['perl', 'score3.pl', 'RecSeqsToAnalyze.fa'])

            os.chdir('../IntronicRiboSnitch')

            # Grab the max ent scores for the sequences

            maxentscores_for_rec = [i.split('\t') for i in pipe_out.strip().split('\n')]

            return maxentscores_for_rec


# This function will return mutations whose mutated sequence has a maxentscore that is comparable to that of the true splice sites
# Remember for calculating max ent score of 5' splice site, you need last three bases of exon and first six bases of intron
# For 3' splice site, you need last 20 bases of intron and first three bases of exon
def returnMutationsCreatingSpliceSitesWithComparableScoresAsTrueSpliceSite(mutseqfile, mutfile, truesplicesitefile,
                                                                        splicesite, min_sites,
                                                                           comparablescore):
    # Grab all the records from mutseqfile store in variable allrecs
    allrecs = []
    for record in SeqIO.parse(open(mutseqfile, 'rU'), 'fasta'):
        allrecs.append(record)

    # Grab all the mutations from the mutfile
    with open(mutfile) as f:
        allmuts = [line.strip().split('\t') for line in f]

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
        rec_mutname = rec.id

        # Since sequence around mutations were stored with a 50 bp window, we only need the 20 base pair window for
        # 5' splice sites
        if splicesite == '5p':
            rec_seq = rec.upper().seq[15:36]
        else:
            rec_seq = rec.upper().seq

        # Find the mutation that matches the rec and if no match is found, break out of loop
        mut_for_rec = [i for i in allmuts if i[3] == rec_mutname]
        #print mut_for_rec
        if len(mut_for_rec) != 1:
            print("No mutation found for record")
            print rec_mutname
            break
        else:
            # If strand is minus, we need to get complement of base and it's mutation, if not just keep base.
            # Use dictionary to grab complement
            if mut_for_rec[0][8].split(';')[3].split('=')[1] == "-":
                mutated_base = base_dict[mut_for_rec[0][5]]
                base_mutated_to = base_dict[mut_for_rec[0][6]]
            else:
                mutated_base = mut_for_rec[0][5]
                base_mutated_to = mut_for_rec[0][6]
            # If the base of concern on the sequence of the record is not the same as the mutated base, break out of loop
            if (splicesite == "5p" and rec_seq[10] != mutated_base) or (
                    splicesite == "3p" and rec_seq[25] != mutated_base):
                print("WT mutation site does not match record")
                print rec_mutname
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
                            seqs_for_rec.append(SeqRecord(seq_i, id=rec_mutname + "_" + str(i + 1) + "_WT"))
                            seqs_for_rec.append(SeqRecord(seq_mutated_i, id=rec_mutname + "_" + str(i + 1) + "_MUT"))
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
                            seqs_for_rec.append(SeqRecord(seq_i, id=rec_mutname + "_" + str(i + 1) + "_WT"))
                            seqs_for_rec.append(SeqRecord(seq_mutated_i, id=rec_mutname + "_" + str(i + 1) + "_MUT"))
                else:
                    if min_sites:
                        for i in range(0, 4):
                            # When i = 0, the mutation site is at last base of the first two bases of the intron,
                            # as i increases, the site will then get closer to the first base of the last two bases
                            # of the exon
                            seq_i = rec_seq[4 + i:27 + i]
                            seq_mutated_i = seq_i[0:21 - i] + base_mutated_to + seq_i[22 - i:]
                            seqs_for_rec.append(SeqRecord(seq_i, id=rec_mutname + "_" + str(i + 1) + "_WT"))
                            seqs_for_rec.append(SeqRecord(seq_mutated_i, id=rec_mutname + "_" + str(i + 1) + "_MUT"))
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
                            seqs_for_rec.append(SeqRecord(seq_i, id=rec_mutname + "_" + str(i + 1) + "_WT"))
                            seqs_for_rec.append(SeqRecord(seq_mutated_i, id=rec_mutname + "_" + str(i + 1) + "_MUT"))

                seqs_for_rec.extend([i for i in truesplicesites if i.id.split('_')[0] == rec_mutname])
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
                    recs_to_consider.append(rec_mutname)
    #if os.getcwd() != '/Users/jayashreekumar/Documents/LaederachLab/Projects/Generating_Statistics_About_Introns':
        #os.chdir('../Projects/Generating_Statistics_About_Introns')
    return recs_to_consider


# This function simply creates new files that contain the mutations that create new splice sites based on MaxEntScan
def filterOutMutationsThatCreateSpliceSites(mutseqfile, mutfile, muts_create5psite_file,muts_create3psite_file,
                                            mutseqfile_new, mutfile_new):
    # Grab all the records from mutseqfile store in variable allrecs
    allrecs = []
    for record in SeqIO.parse(open(mutseqfile, 'rU'), 'fasta'):
        allrecs.append(record)

    # Grab all the mutations from the mutfile
    with open(mutfile) as f:
        allmuts = [line.strip().split('\t') for line in f]

    mutations_createSpliceSites = []
    # Grab the mutations that are predicted to create a stronger splice site than true splice site
    with open(muts_create5psite_file) as f:
        for line in f:
            mutations_createSpliceSites.append(line.strip())
    with open(muts_create3psite_file) as f:
        for line in f:
            mutations_createSpliceSites.append(line.strip())

    # This variable stores fasta records that will be written that don't create splice sites
    recs_tobewritten = [i for i in allrecs if i.id not in mutations_createSpliceSites]

    SeqIO.write(recs_tobewritten, mutseqfile_new, 'fasta')

    # This variable stores bedfile records that will be written that don't create splice sites
    muts_tobewritten = [i for i in allmuts if i[3] not in mutations_createSpliceSites]

    with open(mutfile_new, 'w') as fw:
        for j in muts_tobewritten:
            fw.write('\t'.join(j))
            fw.write('\n')


# This function grab the coordinates for the 2 splice sites surrounding the mutation
# grab sequence last three bases of exon and first six bases of intron and attach it to the mutation for a 5' splice site
# for 3' splice site, grab last 20 bases of intron and first three bases of exon
# Write into two different files for each mutation
# NOTE: sometimes you may get the same splice site being assigned twice to the same mutation if the 5'end or 3' end of
# introns of differing lengths containing the mutation are the same
# It does not make a difference when comparing the true splice site score because they will get the same score since its
# the same site
def assignTrueSpliceSiteClosestToMutationOLDFUNCTION(intersectingfile, strand, splicesite, seqdirectory, writefile):
    # Grab all the mutations in the intersecting file with all column values
    with open(intersectingfile) as f:
        allmuts = [line.strip().split('\t') for line in f]

    # List to store values to write
    muts_truesplicesite_coords = []

    # Grab the mutation ids
    muts_ids = list(set([i[3] for i in allmuts]))

    # Go through each mutation
    for mut in muts_ids:
        # For each mutation grab the introns overlapping that mutation as long as strand that mutation is on is
        # the same as the intron, the mutation is greater than 6 bases away and less than 400 bases from start
        # of the intron.
        # if strand == "plus" and splicesite == "5p":
        #     allintrons_for_mut = [i for i in allmuts if
        #                           i[3] == mut and i[8].split(";")[3][7] == i[14] and int(i[2]) > int(i[10]) + 6 and int(
        #                               i[2]) <= int(i[10]) + 400]
        # elif strand == "minus" and splicesite == "5p":
        #     allintrons_for_mut = [i for i in allmuts if
        #                           i[3] == mut and i[8].split(";")[3][7] == i[14] and int(i[2]) < int(i[11]) - 6 and int(
        #                               i[2]) >= int(i[11]) - 400]
        # elif strand == "plus" and splicesite == "3p":
        #     allintrons_for_mut = [i for i in allmuts if
        #                           i[3] == mut and i[8].split(";")[3][7] == i[14] and int(i[2]) < int(
        #                               i[11]) - 20 and int(i[2]) >= int(i[11]) - 400]
        # else:
        #     allintrons_for_mut = [i for i in allmuts if
        #                           i[3] == mut and i[8].split(";")[3][7] == i[14] and int(i[2]) > int(
        #                               i[10]) + 20 and int(i[2]) <= int(i[10]) + 400]
        allintrons_for_mut = [i for i in allmuts if i[3] == mut]

        # Grab just the unique coordinates for introns.
        coords_allintrons_for_mut = [list(line) for line in set(tuple(i[9:12]) for i in allintrons_for_mut)]

        # For each intron that is valid, grab coordinates of true splice site such that we get three last bases
        # of exon and first six bases of intron
        for i in range(0, len(coords_allintrons_for_mut)):
            if strand == "plus" and splicesite == "5p":
                muts_truesplicesite_coords.append(
                    [mut, str(i + 1), coords_allintrons_for_mut[i][0], str(int(coords_allintrons_for_mut[i][1]) - 2),
                     str(int(coords_allintrons_for_mut[i][1]) + 6)])
            elif strand == "minus" and splicesite == "5p":
                muts_truesplicesite_coords.append(
                    [mut, str(i + 1), coords_allintrons_for_mut[i][0], str(int(coords_allintrons_for_mut[i][2]) - 5),
                     str(int(coords_allintrons_for_mut[i][2]) + 3)])
            elif strand == "plus" and splicesite == "3p":
                muts_truesplicesite_coords.append(
                    [mut, str(i + 1), coords_allintrons_for_mut[i][0], str(int(coords_allintrons_for_mut[i][2]) - 19),
                     str(int(coords_allintrons_for_mut[i][2]) + 3)])
            else:
                muts_truesplicesite_coords.append(
                    [mut, str(i + 1), coords_allintrons_for_mut[i][0], str(int(coords_allintrons_for_mut[i][1]) - 2),
                     str(int(coords_allintrons_for_mut[i][1]) + 20)])
    # Now for those coordinates, grab sequences of true splice site and write to a file
    chrms = range(1, 23)
    chrms.extend(['X', 'Y'])
    # print(chrms)
    chrms_act = ['chr' + str(i) for i in chrms]

    # This variable will contain all true splicing sites
    seqtowrite = []

    for chromosome in chrms_act:
        # print(chromosome)
        mutstograb = [i for i in muts_truesplicesite_coords if i[2] == chromosome]
        # print('Here')
        recsforchrom = []
        for record in SeqIO.parse(open(seqdirectory + '/' + chromosome + '.fa', 'rU'), 'fasta'):
            recsforchrom.append(record)
        for mut in mutstograb:
            if strand == "minus":
                seq = SeqRecord(Seq(str(recsforchrom[0].seq)[int(mut[3]) - 1:int(mut[4])],
                                    IUPAC.unambiguous_dna).reverse_complement(), id=mut[0] + '_' + mut[1],
                                description=mut[0] + ', ' + chromosome + ', start = ' + mut[3] + ', end = ' + mut[4])
            else:
                seq = SeqRecord(Seq(str(recsforchrom[0].seq)[int(mut[3]) - 1:int(mut[4])], IUPAC.unambiguous_dna),
                                id=mut[0] + '_' + mut[1],
                                description=mut[0] + ', ' + chromosome + ', start = ' + mut[3] + ', end = ' + mut[4])

            seqtowrite.append(seq)
            # print('Here3')

    SeqIO.write(seqtowrite, writefile, 'fasta')

def assignTrueSpliceSiteClosestToMutation(intersectingfile, splicesite,seqdirectory, writefile):

    # Grab all the mutations in the intersecting file with all column values
    allmuts = pd.read_csv(intersectingfile,header=None,sep="\t")
    # List to store values to write
    muts_truesplicesite_coords = []

    # Grab the mutation ids
    muts_ids = list(set(allmuts[8]))

    # Go through each mutation
    for mut in muts_ids:
        intron_for_mut = allmuts[allmuts[8]==mut]
        #print mut
        chrom_for_mut = intron_for_mut.iloc[0,0]
        # For each intron that is valid, grab coordinates of true splice site such that we get three last bases
        # of exon and first six bases of intron
        if intron_for_mut.iloc[0,4] == '+':
            if splicesite =='5p':
                muts_truesplicesite_coords.append([mut,chrom_for_mut,intron_for_mut.iloc[0,1]-3,intron_for_mut.iloc[0,1]+6,'+'])
            else:
                muts_truesplicesite_coords.append([mut,chrom_for_mut,intron_for_mut.iloc[0,2]-20,intron_for_mut.iloc[0,2]+3,'+'])
        else:
            if splicesite == '5p':
                muts_truesplicesite_coords.append([mut,chrom_for_mut, intron_for_mut.iloc[0,2]-6, intron_for_mut.iloc[0,2]+3,'-'])
            else:
                muts_truesplicesite_coords.append([mut,chrom_for_mut, intron_for_mut.iloc[0,1]-3 , intron_for_mut.iloc[0,1]+20,'-'])
    # Now for those coordinates, grab sequences of true splice site and write to a file
    chrms = range(1, 23)
    chrms.extend(['X', 'Y'])
    # print(chrms)
    chrms_act = ['chr' + str(i) for i in chrms]

    # This variable will contain all true splicing sites
    seqtowrite = []

    for chromosome in chrms_act:
        print(chromosome)
        mutstograb_5p = [i for i in muts_truesplicesite_coords if i[1] == chromosome]
        # print('Here')
        recsforchrom = []
        for record in SeqIO.parse(open(seqdirectory + '/' + chromosome + '.fa', 'rU'), 'fasta'):
            recsforchrom.append(record)
        for mut in mutstograb_5p:
            if mut[4]=='-':
                seq = SeqRecord(Seq(str(recsforchrom[0].seq)[mut[2]:mut[3]],
                                    IUPAC.unambiguous_dna).reverse_complement().upper(), id=mut[0],
                                description=chromosome + ', start = ' + str(mut[2]+1) + ', end = ' + str(mut[3]) + ', strand=-')
            else:
                seq = SeqRecord(Seq(str(recsforchrom[0].seq)[mut[2]:mut[3]], IUPAC.unambiguous_dna).upper(),
                                id=mut[0],
                                description=chromosome + ', start = ' + str(mut[2]+1) + ', end = ' + str(mut[3]) + ', strand=+')

            seqtowrite.append(seq)
            # print('Here3')

    SeqIO.write(seqtowrite, writefile, 'fasta')
# THis function will run SNPfold on every mutation and write a file that contains the results with the mutation information
def runSNPFold(mutseqfile, mutfile, writefile):
    with open(mutfile) as f:
        mutations = [line.strip().split('\t') for line in f]

    # Create a dictionary that stores the WT and MUT base as values with mutation id as the key
    # If strand is minus, then get complement of bases
    muts_dict = {}
    for i in mutations:
        if i[8].split(';')[3].split('=')[1] == "-":
            muts_dict[i[3]] = [base_dict[i[5]], base_dict[i[6]]]
        else:
            muts_dict[i[3]] = [i[5], i[6]]

    # Grab all the records from mutseqfile store in variable allrecs
    allrecs = []
    for record in SeqIO.parse(open(mutseqfile, 'rU'), 'fasta'):
        allrecs.append(record)

    # Go through every record and run SNPFold on it, the accurate version and store results in variable snpfold_score
    snpfold_score = []
    for rec in allrecs:
        mut_of_interest = muts_dict[rec.id]

        pipe_out = subprocess.check_output(
            ['SNPfold_commandline.py', "-A", str(rec.seq), mut_of_interest[0] + "51" + mut_of_interest[1]])

        snpfold_for_rec = [i.split('\t') for i in pipe_out.strip().split('\n')][1]

        snpfold_score.append([rec.id, snpfold_for_rec])

    with open(writefile, "w") as fw:
        fw.write("MUT" + "\t" + "CC_BPPROB" + "\t" + "CC_BPPROB_PVAL" + "\t" + "CC_SHANNON" + "\t" + "CC_SHANNON_PVAL")
        fw.write("\n")
        for i in snpfold_score:
            fw.write(i[0] + "\t" + "\t".join(i[1][1:]))
            fw.write("\n")

            # return snpfold_score


# THis function will run remuRNA on every mutation and write a file that contains the results with the mutation information
def run_remuRNA(mutseqfile, mutfile, writefile):
    with open(mutfile) as f:
        mutations = [line.strip().split('\t') for line in f]

    # Create a dictionary that stores the WT and MUT base as values with mutation id as the key
    # If strand is minus, then get complement of bases
    muts_dict = {}
    for i in mutations:
        if i[8].split(';')[3].split('=')[1] == "-":
            muts_dict[i[3]] = [base_dict[i[5]], base_dict[i[6]]]
        else:
            muts_dict[i[3]] = [i[5], i[6]]

    # Grab all the records from mutseqfile store in variable allrecs
    allrecs = []
    for record in SeqIO.parse(open(mutseqfile, 'rU'), 'fasta'):
        allrecs.append(record)

    # Go through every record and run remuRNA on it, store results in variable remuRNA_score
    remuRNA_score = []
    for rec in allrecs:
        mut_of_interest = muts_dict[rec.id]

        SeqIO.write(rec, 'RecSeqToAnalyze_remuRNA.fa', 'fasta')

        with open('RecSeqToAnalyze_remuRNA.fa', 'a') as fw:
            fw.write("*" + mut_of_interest[0] + "101" + mut_of_interest[1])

        pipe_out = subprocess.check_output(['/home/mcorley/programs/remuRNA/remuRNA', 'RecSeqToAnalyze_remuRNA.fa'])

        remuRNA_for_rec = [i.split('\t') for i in pipe_out.strip().split('\n')][2]

        remuRNA_score.append([rec.id, remuRNA_for_rec])

    with open(writefile, "w") as fw:
        fw.write("MUT" + "\t" + "MFE(wt)" + "\t" + "MFE(mut)" + "\t" + "dMFE" + "\t" + "H(wt||mu)" + "\t" + "GCratio")
        fw.write("\n")
        for i in remuRNA_score:
            fw.write(i[0] + "\t" + "\t".join(i[1][1:]))
            fw.write("\n")

    return remuRNA_score


# THis function will run RNAsnp on every mutation and write a file that contains the results with the mutation information
def runRNAsnp(mutseqfile, mutfile, writefile):
    with open(mutfile) as f:
        mutations = [line.strip().split('\t') for line in f]

    # Create a dictionary that stores the WT and MUT base as values with mutation id as the key
    # If strand is minus, then get complement of bases
    muts_dict = {}
    for i in mutations:
        if i[8].split(';')[3].split('=')[1] == "-":
            muts_dict[i[3]] = [base_dict[i[5]], base_dict[i[6]]]
        else:
            muts_dict[i[3]] = [i[5], i[6]]

    # Grab all the records from mutseqfile store in variable allrecs
    allrecs = []
    for record in SeqIO.parse(open(mutseqfile, 'rU'), 'fasta'):
        allrecs.append(record)

    # Go through every record and run remuRNA on it, store results in variable remuRNA_score
    RNAsnp_score = []
    for rec in allrecs:
        mut_of_interest = muts_dict[rec.id]

        SeqIO.write(rec, 'RecSeqToAnalyze_RNAsnp.fa', 'fasta')

        with open('SNPToAnalyze_RNAsnp.txt', 'w') as fw:
            fw.write(mut_of_interest[0] + "101" + mut_of_interest[1])

        pipe_out = subprocess.check_output(
            ['/home/mcorley/programs/RNAsnp-1.1/Progs/RNAsnp', "-w", "100", "-l", "10", "-c", "0.0", "-f", 'RecSeqToAnalyze_RNAsnp.fa', "-s",
             "SNPToAnalyze_RNAsnp.txt"])

        RNAsnp_for_rec = [i.split('\t') for i in pipe_out.strip().split('\n')][3]

        RNAsnp_score.append([rec.id, RNAsnp_for_rec])

    with open(writefile, "w") as fw:
        fw.write(
            "MUT" + "\t" + "w" + "\t" + "Slen" + "\t" + "GC" + "\t" + "interval" + "\t" + "dmax" + "\t" + "p-value" + "\t" + "interval" + "\t" + "r_min" + "\t" + "p-value")
        fw.write("\n")
        for i in RNAsnp_score:
            fw.write(i[0] + "\t" + "\t".join(i[1][1:]))
            fw.write("\n")

    return RNAsnp_score


def combineRNAfoldingScoresIntoOneFile(snpfoldfile, remurnafile, rnasnpfile, writefile):
    with open(snpfoldfile) as f:
        lines = [line.strip().split('\t') for line in f][1:]
    snpfolddata = {i[0]: i[1:] for i in lines}

    with open(remurnafile) as f:
        lines = [line.strip().split('\t') for line in f][1:]
    remurnadata = {i[0]: i[1:] for i in lines}

    with open(rnasnpfile) as f:
        lines = [line.strip().split('\t') for line in f][1:]
    rnasnpdata = {i[0]: i[1:] for i in lines}

    with open(writefile, 'w') as fw:
        fw.write(
            "MUT" + "\t" + "SNPfold_CC_BPPROB" + "\t" + "remuRNA_Relative_Entropy" + "\t" + "RNAsnp_dmax" + "\t" + "RNAsnp_dmaxPval" + "\t" + "RNAsnp_rmin" + "\t" + "RNAsnp_rminPval" + "\n")
        for key in snpfolddata:
            fw.write(key + "\t" + snpfolddata[key][0] + "\t" + remurnadata[key][3] + "\t" + rnasnpdata[key][4] + "\t" +
                     rnasnpdata[key][5] + "\t" + rnasnpdata[key][7] + "\t" + rnasnpdata[key][8] + "\n")


def findMutationsNearbyCurrentMut(mutfile, allmutsfile, window, writefile):
    with open(mutfile) as f:
        muts = [line.strip().split('\t') for line in f]

    with open(allmutsfile) as f:
        allmuts = [line.strip().split('\t') for line in f]
    with open(writefile, 'w') as fw:
        for m in muts:
            name_m = m[3]
            mut_coord = int(m[2])
            mut_chr = m[0]
            strand = m[8].split(";")[3][7]
            mutsnearby = [i for i in allmuts if i[0] == mut_chr and i[8].split(";")[3][7] == strand and int(
                i[2]) > mut_coord - window and int(i[2]) < mut_coord + window]

            fw.write(name_m + "\t" + str(len(mutsnearby)) + "\n")


def plotSNPFold_BPprobability(allbpprob, window):
    with open(allbpprob) as f:
        all = [line.strip().split('\t') for line in f][1:]

    WT_score = [float(i[1]) for i in all] #if int(i[0]) >= 115 - window and int(i[0]) <= 115 + window]
    MUT_score = [float(i[2]) for i in all]#if int(i[0]) >= 115 - window and int(i[0]) <= 115 + window]
    print len(WT_score)
    print len(MUT_score)
    print len(range(0, 2 * window + 1))
    plt.plot(range(0, 2 * window + 1), MUT_score, color="red", label="MUT")
    plt.plot(range(0, 2 * window + 1), WT_score, color="black", label="WT")
    plt.xlabel("Nucleotide Position")
    #y = range(-window, window)
    plt.xlim(0,2*window+1)
    x = range(0, 2 * window + 1,20)
    plt.xticks(x,range(0, 2 * window + 1,20))

    plt.ylabel("Base-pair probability")
    plt.axvline(x=window, ymin=0, ymax=1, color="green")
    plt.legend(loc='upper left')
    plt.show()


def calculateAndPlotDistributionOfMutsAroundAndWithinISREs(mutationsWithMotifInfo,motifFile):

    with open(motifFile) as f:
        ISREs = [line.strip().split('\t') for line in f]

    pos_for_all_ISREs = []
    for m in ISREs:
        pos_for_m = []
        for f in mutationsWithMotifInfo:
            for i in f[2]:
                if i[0] == m:
                    for j in i[1]:
                        pos_for_m.append(j)
        pos_for_all_ISREs.append([m, pos_for_m])

    combinedpos_ISREs = []
    for i in pos_for_all_ISREs:
        combinedpos_ISREs.extend(i[1])

    freq_pos_ISREs = Counter(combinedpos_ISREs)
    freq_pos_ISREs_list = [[k, v] for k, v in freq_pos_ISREs.iteritems()]

    freq_pos_ISRE_list_sorted = sorted(freq_pos_ISREs_list, key=itemgetter(0))
    plt.bar([i[0] for i in freq_pos_ISRE_list_sorted], [i[1] for i in freq_pos_ISRE_list_sorted], align='center')
    plt.xticks([i[0] for i in freq_pos_ISRE_list_sorted])



def calculateIfMFEChangesWithMutation(mutationsOfInterest,mutationfile,writefile):

    #allmutseq = dict()

    #for record in SeqIO.parse(open(mutationseqfile, 'rU'), 'fasta'):
        #allmutseq[record.id] = record
    with open(mutationfile) as f:
        allmutations = [line.strip().split('\t') for line in f]
    mut_dict = dict()
    for mutation in allmutations:
        mut_dict[mutation[3]] = [mutation[5],mutation[6],mutation[8].split(';')[3].split('=')[1]]
    with open(writefile,'w') as fw:
        for mut in mutationsOfInterest[0:2]:
            seq = str(mut[1])
            mutid = mut[0]
            mut_info = mut_dict[mutid]
            for i in mut[2]:
                motif = i[0][0]
                for j in i[1]:
                    # first index is 25-j+10, second index is 25-j+5+10+1
                    WTseq = SeqRecord(Seq(seq[15-j:41-j],IUPAC.unambiguous_dna),id= mutid+'_'+motif+'_WT')
                    SeqIO.write([WTseq], "FileToSendToFold.fa", 'fasta')
                    pipe_out = subprocess.check_output(['Fold', 'FileToSendToFold.fa', 'CTfile.ct','--loop', '30','--maximum','20','--percent', '10','--temperature','310.15','--window','3'])
                    with open('CTfile.ct') as f:
                        lines = [line.strip() for line in f]
                    WT_energy = lines[0].split(' ')[4]
                    if mut_info[2] == "-" and base_dict[mut_info[0]] == seq[25]:
                        MUTseq = SeqRecord(Seq(seq[15-j:25] + base_dict[mut_info[1]] + seq[26:41-j],IUPAC.unambiguous_dna),id= mutid+'_'+motif+'_MUT')
                        SeqIO.write([MUTseq], "FileToSendToFold.fa", 'fasta')
                        pipe_out = subprocess.check_output(
                            ['Fold', 'FileToSendToFold.fa', 'CTfile.ct', '--loop', '30', '--maximum', '20', '--percent',
                             '10', '--temperature', '310.15', '--window', '3'])
                        with open('CTfile.ct') as f:
                            lines = [line.strip() for line in f]

                        MUT_energy = lines[0].split(' ')[4]
                    elif mut_info[2] == "+" and mut_info[0] == seq[25]:
                        MUTseq = SeqRecord(Seq(seq[15-j:25] + mut_info[1] + seq[26:41-j],IUPAC.unambiguous_dna),id= mutid+'_'+motif+'_MUT')
                        SeqIO.write([MUTseq], "FileToSendToFold.fa", 'fasta')
                        pipe_out = subprocess.check_output(
                            ['Fold', 'FileToSendToFold.fa', 'CTfile.ct', '--loop', '30', '--maximum', '20', '--percent',
                             '10', '--temperature', '310.15', '--window', '3'])
                        with open('CTfile.ct') as f:
                            lines = [line.strip() for line in f]

                        MUT_energy = lines[0].split(' ')[4]
                    else:
                        print "Base does not match for " + mutid
                        break

                    fw.write(mutid + '\t'+ motif + '\t' + str(j)+'\t'+str(WT_energy)+'\t'+str(MUT_energy))

# This function will create a fasta file that includes the mutated base within the Wt sequence
def createFastaFilewithMutationsIncluded(mutbedfile,wtfastafile,strand,mutfastawritefile):

    # Open the mutation file which includes the wt base and the mutated base and create a dictionary that uses
    # the mutation id as the key and the wt base and mutated base as the values
    # Need to do reverse complement if the mutation is on the minus strand
    with open(mutbedfile) as f:
        allmutations = [line.strip().split('\t') for line in f]
    mut_dict = dict()
    for mutation in allmutations:
        if strand == 'plus':
            mut_dict[mutation[3]] = [mutation[5],mutation[6]]
        else:
            mut_dict[mutation[3]] = [base_dict[mutation[5]], base_dict[mutation[6]]]

    # Collect all the new records which includes the mutated base within the wt sequence
    # So go through every record in the fasta file, and for that record id, find the wt base and mutated base and
    # in the fasta sequence replace the WT base with the MUT base.
    new_records = []
    for record in SeqIO.parse(open(wtfastafile, 'rU'), 'fasta'):
        seq_for_record = str(record.seq).upper()
        mut_for_record = mut_dict[record.id]
        if mut_for_record[0] == seq_for_record[500]:
            seq = SeqRecord(
                Seq(seq_for_record[0:500]+mut_for_record[1]+seq_for_record[501:], IUPAC.unambiguous_dna),
                id=record.id,
                description=record.description)
            new_records.append(seq)
        else:
            print "Base does not match for " + record.id

    # Write new records to a new fasta file
    SeqIO.write(new_records, mutfastawritefile, 'fasta')


def getBasePairingBPsForWTandMUT(mutdata,mutseq):
    # Create a dictionary that stores the WT and MUT base as values with mutation id as the key
    # If strand is minus, then get complement of bases

    if mutdata[8].split(';')[3].split('=')[1] == "-":
        mut_of_interest = [base_dict[mutdata[5]], base_dict[mutdata[6]]]
    else:
        mut_of_interest = [mutdata[5], mutdata[6]]

    #pipe_out = subprocess.check_output(
    #        ['SNPfold_commandline.py', "-A", "-n",mutdata[3],"-m","bpProbs", str(mutseq), mut_of_interest[0] + "51" + mut_of_interest[1]])

    bp_probs_all = pd.read_csv("output/"+mutdata[3]+"/input_mutations_bpProbs_perNt.txt",sep="\t",header=None,index_col=0)
    WT_bp_probs = bp_probs_all.loc["WT"]
    if mut_of_interest[0] == 'T':
        mut_of_interest[0] = 'U'
    if mut_of_interest[1] == 'T':
        mut_of_interest[1] = 'U'
    MUT_bp_probs = bp_probs_all.loc[mut_of_interest[0] + "51" + mut_of_interest[1]]

    return [WT_bp_probs,MUT_bp_probs]

def mutationsThatAlterBPOfMotifs(datafile,mutfile,mutseqfile,bpdiff,writefile):

    databeforegroup = pd.read_csv(datafile,header=None,sep="\t")
    data = databeforegroup.groupby(8)

    with open(mutfile) as f:
        mutations = [line.strip() for line in f]

    recs_dict = {}
    for record in SeqIO.parse(open(mutseqfile, 'rU'), 'fasta'):
        recs_dict[record.id] = record

    with open(writefile,'w') as fw:
        for mut in mutations:
            print mut
            data_for_mut = data.get_group(mut)
            mut_data = list(data_for_mut.iloc[0,range(5,14)])
            mut_seq = str(recs_dict[mut].seq)
            bp_wtandmut = getBasePairingBPsForWTandMUT(mut_data,mut_seq)
            data_start_sites = list(data_for_mut[6]-(data_for_mut[1]+100))
            fw.write(mut + '\t' + str(len(data_start_sites)))
            for i in data_start_sites:
                if i < 0:
                    print 'Here'
                    bp_diff = abs(bp_wtandmut[0][50+abs(i):abs(i)+56] - bp_wtandmut[1][50+abs(i):abs(i)+56])
                else:
                    bp_diff = abs(bp_wtandmut[0][50-i:56-i] - bp_wtandmut[1][50-i:56-i])
                #bp_prob_same = sum(bp_diff == 0)
                fw.write('\t' + str(sum(bp_diff>bpdiff)))
            fw.write('\n')