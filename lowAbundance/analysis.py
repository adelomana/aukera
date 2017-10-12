###
### this script evaluates if low abundant transcripts whose proteins upregulate, have specific footprint patterns
###

import sys
import matplotlib,matplotlib.pyplot
matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})

def proteomicsReader():

    data={} # data[lysate/rbf][rep1/rep2/rep3][tp2vs1/tp3vs1/tp4vs1][geneName]=log2FC
    significance={} # same as data
    
    allFiles=os.listdir(proteomicsDataFolder)
    csvFiles=[element for element in allFiles if '.csv' in element and '._' not in element]

    for csvFile in csvFiles:
        path=proteomicsDataFolder+csvFile

        brokenName=csvFile.split('.')
        condition=brokenName[0]
        replicate=brokenName[1]

        if condition not in data.keys():
            data[condition]={}; significance[condition]={}
        if replicate not in data[condition].keys():
            data[condition][replicate]={}; significance[condition][replicate]={}
            
        timepoints=['tp2vs1','tp3vs1','tp4vs1']
        for timepoint in timepoints:
            if timepoint not in data[condition][replicate].keys():
                data[condition][replicate][timepoint]={}; significance[condition][replicate][timepoint]={}

        with open(path,'r') as f:
            next(f)
            for line in f:
                vector=line.split(',')
                geneName=vector[0]

                a=float(vector[2])
                b=float(vector[6])
                c=float(vector[10])

                d=float(vector[2+2])
                e=float(vector[6+2])
                f=float(vector[10+2])

                data[condition][replicate]['tp2vs1'][geneName]=a
                data[condition][replicate]['tp3vs1'][geneName]=b
                data[condition][replicate]['tp4vs1'][geneName]=c

                significance[condition][replicate]['tp2vs1'][geneName]=d
                significance[condition][replicate]['tp3vs1'][geneName]=e
                significance[condition][replicate]['tp4vs1'][geneName]=f
            
    return data,significance

###
### MAIN
###

# 0. user defined variables


# 0.1. paths
proteomicsDataFolder='/Volumes/omics4tb/alomana/projects/TLR/data/proteomics/all/'

# 0.2. variables
timepoints=[14.3,21.5,28.8,40.8] 
sortedReplicateLabels=['br1','br2','br3']
sortedRelativeTimePointLabels=['tp2vs1','tp3vs1','tp4vs1']

# 1. reading data
# 1.1. reading protein data
log2proteome,proteomeSignificance=proteomicsReader()

# define all up-regulated proteins

# intersect with the 25% lower abundant transcripts

# remove all candidate genes that are transcriptionally upregulated, in other words, the pt change can be explained by transcript changes

# inspect the footprint load in those. 
