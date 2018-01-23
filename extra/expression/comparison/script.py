### this script compares expression of kallisto/cufflinks/HTseq-count using rank correlations

#### right now it's more about comparing my counts with Arjuns counts

import sys,os

def readerADL():

    '''
    This function reads ADL counts.
    '''

    countsADL={}
    conditions=[]; geneNames=[]

    dataFiles=os.listdir(countsADLdir)
    for dataFile in dataFiles:
        fileName=countsADLdir+dataFile
        flag=dataFile.split('.txt')[0]

        countsADL[flag]={}

        if flag not in conditions:
            conditions.append(flag)

        with open(fileName,'r') as f:
            for line in f:
                vector=line.split('\t')
                geneName=vector[0]
                value=float(vector[1])
                countsADL[flag][geneName]=value

                if geneName not in geneNames:
                    geneNames.append(geneName)

    conditions.sort()
    geneNames.sort()

    return countsADL,conditions,geneNames

def readerAR():

    '''
    This function reads AR counts.
    '''

    countsAR={}

    dataFiles=os.listdir(countsARdir)

    for dataFile in dataFiles:
        fileName=countsARdir+dataFile

        if 'mrna' in dataFile:
            method='trna'
        elif 'rbf' in dataFile:
            method='rbf'
        else:
            print('problem a')
            sys.exit()

        with open(fileName,'r') as f:
            firstLine=f.readline()
            header=firstLine.split(',')
            
            for line in f:
                vector=line.split(',')

                # geneName
                geneName=vector[0]

                # define timepoint, replicate and method
                for i in range(len(header)):
                    print(i,header[i])

                    
                sys.exit()
    

    return countsAR


###
### MAIN
###

# 0. user defined variables
countsADLdir='/Volumes/omics4tb/alomana/projects/TLR/data/counts/'
countsARdir='/Volumes/omics4tb/alomana/projects/TLR/data/countsArjun/'

# 1. read the data
countsADL,conditions,geneNames=readerADL()
countsAR=readerAR()

# 2. make the comparisons
