import pandas as pd
import subprocess
# This function will run SNPfold on every mutation and write a file that contains the results with the mutation information
# It requires a file that contains sequences around mutation, a file containing info of the WT base and MUT base
def runSNPFold(mutseqfile, mutfile, baseOfinterest,writefile):
    
    # Read the file containing sequences around mutation
    mutseqfile = pd.read_table(mutseqfile,sep="\t",header=None)
    # Create a dictionary for those sequences using the mutation id as key
    mutseq_dict = pd.Series(mutseqfile[1].values,index=mutseqfile[0].values).to_dict()
    
    # Read the file containing information about WT base and MUT base
    mutfile = pd.read_table(mutfile,sep="\t",header=None).iloc[:,range(5,13)].drop_duplicates()
    print mutfile.head()
    # Create a dictionary for strand, WT base and MUT base using the mutation id as key
    mut_dict = {row[8]:[row[9],row[10],row[11]] for index, row in mutfile.iterrows()} 
    #pd.Series(mutfile.iloc[:,[4,5,6]].values,index=mutfile[8].values).to_dict()

    # Dictionary for complementary base dictionary
    base_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    
    # Go through every mutation and run SNPFold on it, the accurate version and write results to file
    with open(writefile, "w") as fw:
        # First write header
        fw.write("MUT" + "\t" + "CC_BPPROB" + "\t" + "CC_BPPROB_PVAL" + "\t" + "CC_SHANNON" + "\t" + "CC_SHANNON_PVAL")
        fw.write("\n")
        # Go through every mutation
        for mut in mut_dict.keys():
            # Get the sequence and bases for the mutation
            seq4mut = mutseq_dict[mut]
            bases4mut = mut_dict[mut]
            # If the mutation is on the negative strand, get the complement of the WT and MUT base
            # Else just get the WT and MUT base
            if bases4mut[0]=="-":
                WTbase=base_dict[bases4mut[1]]
                MUTbase=base_dict[bases4mut[2]]
            else:
                WTbase=bases4mut[1]
                MUTbase=bases4mut[2]
            # Call SNP fold and write result into file    
            try:
	        pipe_out = subprocess.check_output(["SNPfold_commandline.py","-m","all","-n",mut,"-A", seq4mut, WTbase + baseOfinterest + MUTbase])
            	snpfold_for_mut = [i.split('\t') for i in pipe_out.strip().split('\n')][1]
            	fw.write(mut + "\t" + "\t".join(snpfold_for_mut[1:]))
            	fw.write("\n")
	    except:
		print mut

import sys

numBasesToFold = int(sys.argv[1])
whichDataSet = "CommonSNPs"

runSNPFold("../temp/"+whichDataSet+"/"+whichDataSet+"_Flanking_100bpISREs_"+str(numBasesToFold)+"bpwindow-Sequences.txt","../temp/"+whichDataSet+"/"+whichDataSet+"_Flanking_100bpISREs.txt",str((numBasesToFold/2)+1),"../temp/"+whichDataSet+"/"+whichDataSet+"_Flanking_100bpISREs_"+str(numBasesToFold)+"bpwindow_SNPFoldResults.txt")
