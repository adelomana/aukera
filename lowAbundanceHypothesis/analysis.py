###
### this script evaluates if low abundant transcripts whose proteins upregulate, have specific footprint patterns
###

import sys,os,numpy
import matplotlib,matplotlib.pyplot
matplotlib.rcParams.update({'font.size':20,'font.family':'Arial','xtick.labelsize':16,'ytick.labelsize':16})

def consistencyNameChecker():

    '''
    This function determines the intersect of protein and transcript names.
    '''

    transcriptomeNames=[]
    for fraction in rnaExpression.keys():
        for replicate in rnaExpression[fraction].keys():
            for timepoint in rnaExpression[fraction][replicate].keys():
                for name in rnaExpression[fraction][replicate][timepoint].keys():
                    if name not in transcriptomeNames:
                        transcriptomeNames.append(name)

    proteomeNames=[]
    for fraction in log2proteome.keys():
        for replicate in log2proteome[fraction].keys():
            for timepoint in log2proteome[fraction][replicate].keys():
                for name in log2proteome[fraction][replicate][timepoint].keys():
                    if name not in proteomeNames:
                        proteomeNames.append(name)

    print('\t found expression quantification for {} transcripts.'.format(len(transcriptomeNames)))
    print('\t found expression quantification for {} proteins.'.format(len(proteomeNames)))

    # define which proteins do not have transcript equivalent
    consistentNames=[]
    inconsistentNames=[]
    for ptName in proteomeNames:
        if ptName in transcriptomeNames:
            consistentNames.append(ptName)
        else:
            inconsistentNames.append(ptName)
    consistentNames.sort()
    print('\t found transcriptome info for {} proteins.'.format(len(consistentNames)))
    print('\t lost {} proteins for annotation discrepancies:'.format(len(inconsistentNames)))
    print('\t\t {}'.format(','.join(inconsistentNames)))
    print()

    return consistentNames

def proteomicsReader():

    '''
    This function reads proteome data.
    '''

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

def transcriptomeRelativeConverter():

    '''
    this function assumes consistent gene names for all conditions
    '''

    log2transcriptome={}
    
    for fraction in rnaExpression.keys():
        for replicate in rnaExpression[fraction].keys():
            for timepoint in rnaExpression[fraction][replicate].keys():
                if timepoint != 'tp.1':
                    for name in rnaExpression[fraction][replicate][timepoint].keys():
                        a=rnaExpression[fraction][replicate][timepoint][name]
                        b=rnaExpression[fraction][replicate]['tp.1'][name]

                        fc=(a+1)/(b+1)
                        log2FC=numpy.log2(fc)

                        stripped=timepoint.replace('.','')
                        relativeTimePoint=stripped+'vs1'

                        if fraction not in log2transcriptome.keys():
                            log2transcriptome[fraction]={}
                        if replicate not in log2transcriptome[fraction].keys():
                            log2transcriptome[fraction][replicate]={}
                        if relativeTimePoint not in log2transcriptome[fraction][replicate].keys():
                            log2transcriptome[fraction][replicate][relativeTimePoint]={}

                        log2transcriptome[fraction][replicate][relativeTimePoint][name]=log2FC

    return log2transcriptome

def transcriptomicsReader():

    '''
    this function reads transcriptomics data as in
    transcriptomics[trna/rbf][replicate][timepoint][gene]
    '''

    data={}
    
    with open(transcriptomicsDataFile,'r') as f:
        header=f.readline()
        labels=header.split('\t')[1:-1]
        for label in labels:
            crumbles=label.split('.')

            fraction=crumbles[0]
            replicate='br'+crumbles[2]
            timepoint='tp.'+crumbles[4]

            if fraction not in data.keys():
                data[fraction]={}
            if replicate not in data[fraction].keys():
                data[fraction][replicate]={}
            if timepoint not in data[fraction][replicate].keys():
                data[fraction][replicate][timepoint]={}
            
        for line in f:
            vector=line.split('\t')[:-1]
            values=[float(element) for element in vector[1:]]
            geneName=vector[0].replace('_','')
            for i in range(len(values)):
                crumbles=labels[i].split('.')
                fraction=crumbles[0]
                replicate='br'+crumbles[2]
                timepoint='tp.'+crumbles[4]                
                data[fraction][replicate][timepoint][geneName]=values[i]
    
    return data

def upregulatedProteinsFinder():

    '''
    This function defines up-regulated proteins. 
    Up-regulated proteins are defined as:
    1) have significant p-value for all 3 biological replicates
    2) log2 FC > 1 
    3) maybe,  ***check*** they are all within 1 std
    '''

    upPt=[]
    condition='lysate'
    replicates=proteomeSignificance[condition].keys()
    timepoint='tp4vs1'

    x=[]; y=[]; xa=[]; ya=[]

    noLastTime=[]
    
    for geneName in consistentNames:
        foldChanges=[]; pvalues=[]
        for replicate in replicates:
            fc=None; pval=None
            try:
                fc=log2proteome[condition][replicate][timepoint][geneName]
                pval=proteomeSignificance[condition][replicate][timepoint][geneName]
            except:
                pass
            if fc != None:
                foldChanges.append(fc)
                pvalues.append(pval)
        if foldChanges != []:
            x.append(numpy.median(foldChanges))
            y.append(numpy.log10(numpy.median(pvalues)))
            if min(foldChanges) > 1 and max(pvalues) < 0.05:
                upPt.append(geneName)
                xa.append(numpy.median(foldChanges))
                ya.append(numpy.log10(numpy.median(pvalues)))
        else:
            noLastTime.append(geneName)
            
    print('missing {} proteins because of no info on last time point...'.format(len(noLastTime)))
    print(noLastTime[:10])
    print('')

    ### 3rd point was also searched for noLastTime proteins--no data available either

    print('{} proteins found up-regulated'.format(len(upPt)))
    print(upPt[:10])
    print('')

    # computing histograms
    matplotlib.pyplot.hist(x,bins=10,range=(-5,5),label='all')
    matplotlib.pyplot.hist(xa,bins=10,range=(-5,5),label='significant')
    matplotlib.pyplot.xlabel('log$_2$ FC')
    matplotlib.pyplot.ylabel('# of proteins')
    matplotlib.pyplot.legend()
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig('histo.x.pdf')
    matplotlib.pyplot.clf()

    matplotlib.pyplot.hist(y,bins=7,range=(-7,0),label='all')
    matplotlib.pyplot.hist(ya,bins=7,range=(-7,0),label='up-regulated')
    matplotlib.pyplot.xlabel('log$_{10}$ p-value')
    matplotlib.pyplot.ylabel('# of proteins')
    matplotlib.pyplot.legend()
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig('histo.y.pdf')
    matplotlib.pyplot.clf()

    return upPt

###
### MAIN
###

# 0. user defined variables

# 0.1. paths
transcriptomicsDataFile='/Volumes/omics4tb/alomana/projects/TLR/data/expression/expressionMatrix.kallisto.txt'
proteomicsDataFolder='/Volumes/omics4tb/alomana/projects/TLR/data/proteomics/all/'

# 0.2. variables
timepoints=[14.3,21.5,28.8,40.8] 
sortedReplicateLabels=['br1','br2','br3']
sortedRelativeTimePointLabels=['tp2vs1','tp3vs1','tp4vs1']

# 1. reading data

# 1.1. reading mRNA data
rnaExpression=transcriptomicsReader()
log2transcriptome=transcriptomeRelativeConverter()

# 1.2. reading protein data
log2proteome,proteomeSignificance=proteomicsReader()

# 1.3. checking consistency of transcript and protein names
consistentNames=consistencyNameChecker()

# 2. define up-regulated proteins
upPt=upregulatedProteinsFinder()

# 3. define expression for up-regulated proteins at initial timepointintersect with the 25% lower abundant transcripts
upPt_rna=[]
all_rna=[]
lowCandidate=[]

fraction='trna'
timepoint='tp.1'
for geneName in consistentNames:
    values=[]
    for replicate in sortedReplicateLabels:
        expressionValue=rnaExpression[fraction][replicate][timepoint][geneName]
        values.append(numpy.log2(expressionValue+1))
    average=numpy.mean(values)
    
    all_rna.append(average)
    if geneName in upPt:
        upPt_rna.append(average)
        if average < 5:
            lowCandidate.append(geneName)

print(min(all_rna),max(all_rna))
# computing histograms
matplotlib.pyplot.hist(all_rna,bins=15,range=(0,15),label='all')
matplotlib.pyplot.hist(upPt_rna,bins=15,range=(0,15),label='up-reg. pt.')
matplotlib.pyplot.xlabel('log$_2$ (TPM+1)')
matplotlib.pyplot.ylabel('# of genes')
matplotlib.pyplot.legend()
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig('histo.expression.pdf')
matplotlib.pyplot.clf()

# show the trends of pt and rna for all low expressed genes [<5 log2(TPM+1)]
print(len(lowCandidate))
for geneName in lowCandidate:
    print(geneName)
    p=[]; r=[] # get protein and RNA trajectories
    for replicate in sortedReplicateLabels:
        pp=[]; rr=[]
        for timepoint in sortedRelativeTimePointLabels:
            if geneName in log2proteome['lysate'][replicate][timepoint].keys():
                value=log2proteome['lysate'][replicate][timepoint][geneName]
            else:
                value=float('nan')
            pp.append(value)
            if geneName in log2transcriptome['trna'][replicate][timepoint].keys():
                value=log2transcriptome['trna'][replicate][timepoint][geneName]
            else:
                value=float('nan')
            rr.append(value)
        print(replicate,'rna',rr,'pt',pp)
        # adding replicates
        p.append(pp)
        r.append(rr)

    # plotting data
    print('')
    R=numpy.transpose(numpy.array(r))
    P=numpy.transpose(numpy.array(p))
    matplotlib.pyplot.plot(timepoints[1:],R,'o',color='C0')
    matplotlib.pyplot.plot(timepoints[1:],P,'o',color='C1')

    # plotting the estimated linear average
    averageR=numpy.nanmean(R,axis=1)
    averageRref=numpy.insert(averageR,0,0)
    averageP=numpy.nanmean(P,axis=1)
    averagePref=numpy.insert(averageP,0,0)
    print('averageR',averageRref)
    print('averageP',averagePref)
    matplotlib.pyplot.plot(timepoints,averageRref,'-',color='C0',lw=2,label='mRNA')
    matplotlib.pyplot.plot(timepoints,averagePref,'-',color='C1',lw=2,label='pt')
    
    # closing figure
    matplotlib.pyplot.legend()
    matplotlib.pyplot.xlabel('Time, hours')
    matplotlib.pyplot.ylabel('log$_2$ FC')
    matplotlib.pyplot.xlim([-1,42])
    matplotlib.pyplot.axvline(timepoints[0],ls='--',color='grey')
    matplotlib.pyplot.title(geneName)
    fileName='figures/trajectory.{}.pdf'.format(geneName)
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(fileName)
    matplotlib.pyplot.clf()

    #sys.exit()
            

# remove all candidate genes that are transcriptionally upregulated, in other words, the pt change can be explained by transcript changes
# how much can I explain pt by footprints or mRNA changes.

# inspect the footprint load in those. 
