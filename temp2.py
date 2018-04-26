# Import relevant packages
import pandas as pd

def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

from operator import itemgetter
def getMAF_Index(rowofinterest):
    x =rowofinterest[8].split(';')
    y = [i for i in x if i is not None and "CAF=" in i]
    z = y[0].split('=')[1].split(',')
    z_p = [float(i) if isfloat(i) else 1 for i in z[1:]]
    min_index = min(enumerate(z_p), key=itemgetter(1))[0]
    return min_index

def getMAF(rowofinterest):
    x =rowofinterest[8].split(';')
    y = [i for i in x if i is not None and "CAF=" in i]
    z = y[0].split('=')[1].split(',')
    z_p = [float(i) for i in z[1:] if isfloat(i)]
    return min(z_p) 

def getRowRelevantForGroup(group):
    i = group["maf_index"].values[0]
    try:
        return group.iloc[i,:]
    except:
        print group

def getNumOfAlleles(rowofinterest):
    x =rowofinterest[8].split(';')
    y = [i for i in x if i is not None and "CAF=" in i]
    z = y[0].split('=')[1].split(',')
    return len(z)

def getGene(rowofinterest):
    x =rowofinterest[8].split(';')
    y = [i for i in x if i is not None and "GENEINFO=" in i]
    return y[0].split('=')[1]

def extractValidInfoCommonSNPs(snpfile,writefile):
    chrom = range(1,23)
    chrom.extend(['X','Y'])
    commonSNPs = pd.read_csv(snpfile,sep="\t",header=None)
    commonSNPs_single = commonSNPs[(commonSNPs[5].isin(['A','G','T','C','N']))&(commonSNPs[6].isin(['A','G','T','C','N']))]
    for chromosome in chrom:
        print chromosome
        # Get chromosomes that are for each chromosome
        commonSNPs_single_chrom = commonSNPs_single[commonSNPs_single[0]==chromosome]
        # Get number of alleles for each SNP
        numAllelesForChrom = commonSNPs_single_chrom.apply(getNumOfAlleles,1)
        # Get SNPs that have 2 alleles, one major and one minor allele
        dataForchrom_2alleles = commonSNPs_single_chrom[numAllelesForChrom == 2]
        # Get the minor allele frequency for SNPs that have 2 alleles
        mafForchrom_2alleles = dataForchrom_2alleles.apply(getMAF,1)
        # Get all SNPs that have greater than 2 alleles 
        not2alleles = commonSNPs_single_chrom[numAllelesForChrom > 2]
        # Get the index of the allele that has the minor allele frequency, this will 0,1,2 or 3
        index_MAF_not2alleles = not2alleles.apply(getMAF_Index,1)
        # Assign the index to the table that is for SNPs with more than 2 alleles
        not2alleles = not2alleles.assign(maf_index=index_MAF_not2alleles)
        # Group SNPs by SNP id and then get the row that has the index for the MAF
        dataForchrom_not2alleles = not2alleles.groupby([3]).apply(getRowRelevantForGroup)
        # get the minor allele frequency for that allele selected 
        maf_not2alleles = dataForchrom_not2alleles.apply(getMAF,1)
        dataForchrom = {"chrom":pd.concat([dataForchrom_2alleles[0],dataForchrom_not2alleles[0]]),
                        "start":pd.concat([dataForchrom_2alleles[1],dataForchrom_not2alleles[1]]),
                        "end":pd.concat([dataForchrom_2alleles[2],dataForchrom_not2alleles[2]]),
                        "snpID":pd.concat([dataForchrom_2alleles[3],dataForchrom_not2alleles[3]]),
                        "major":pd.concat([dataForchrom_2alleles[5],dataForchrom_not2alleles[5]]),
                        "minor":pd.concat([dataForchrom_2alleles[6],dataForchrom_not2alleles[6]]),
                        "MAF":pd.concat([mafForchrom_2alleles,maf_not2alleles])}
        #dfForchrom = pd.DataFrame(dataForchrom,columns=["chrom","start","end","snpID","major","minor","gene","MAF"])
        dfForchrom = pd.DataFrame(dataForchrom,columns=["chrom","start","end","snpID","major","minor","MAF"])
        dfForchrom.to_csv(writefile,sep="\t",header=False,index=False,mode="a")
        #dfForchrom_MAFvalid = dfForchrom[(dfForchrom["MAF"]>=minMAFval)&(dfForchrom["MAF"]<=(1-minMAFval))]
        #dfForchrom_MAFvalid['chrom'] = 'chr' + dfForchrom_MAFvalid['chrom'].astype(str)
        #dfForchrom_MAFvalid.to_csv(writefile,sep="\t",header=False,index=False,mode="a")

extractValidInfoCommonSNPs("../processed_data/common_all.bed","../processed_data/commonSNPs_all_processed.bed")
