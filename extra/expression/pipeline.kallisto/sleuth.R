###
### this script defines differentially expressed transcripts using Sleuth
###

library('sleuth')

# 0. user defined variables
resultsDir='/Volumes/omics4tb/alomana/projects/TLR/data/kallisto1e2'
metadataFileRBF='metadata.sleuth.csv'

setwd("~/github/aukera/extra/expression/pipeline.kallisto")

# 1. reading data
sampleIDs=dir(file.path(resultsDir))
dirs=file.path(resultsDir,sampleIDs)

# 2. reading metadata
s2c=read.table(metadataFileRBF,header=TRUE,stringsAsFactors=FALSE,sep=",")
#s2c=dplyr::select(s2c, sample=sampleIDs,condition)
