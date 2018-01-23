###
### this script defines differentially expressed transcripts using DESeq2
###

library('DESeq2')
library('stringr')

# 0. user defined variables
setwd("~/github/aukera/extra/expression/pipeline.star.htseq-count.deseq2")
#resultsDir='/Volumes/omics4tb/alomana/projects/TLR/data/kallisto1e3'
countsDir="/Volumes/omics4tb/alomana/projects/TLR/data/counts"

# 1. handle data reading
sampleIDs=dir(file.path(countsDir))

# 1.1. working with RBF
sampleSubset=sampleIDs[str_detect(sampleIDs,"rbf")]
t1=sampleSubset[str_detect(sampleSubset,"tp.1")]
t4=sampleSubset[str_detect(sampleSubset,"tp.4")]
samplesRBF=c(t1,t4)
filesRBF=file.path(countsDir,samplesRBF)

timeCondition=sub("(.*tp).*","\\12",samplesRBF)

# 2. formatting variables
sampleTable=data.frame(sampleName=filesRBF,fileName=filesRBF,condition=samplesRBF)
