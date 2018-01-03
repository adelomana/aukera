###
### this script defines differentially expressed transcripts using Sleuth
###

library('sleuth')
library('stringr')

# 0. user defined variables
setwd("~/github/aukera/extra/expression/pipeline.kallisto")
resultsDir='/Volumes/omics4tb/alomana/projects/TLR/data/kallisto1e2'
metadataFileRBF='sleuth.metadata.rbf.csv'
metadataFileRNA='sleuth.metadata.rna.csv'

# 1. handle data reading
sampleIDs=dir(file.path(resultsDir))
kallistoDirs=file.path(resultsDir,sampleIDs)

sampleIDsRBF=sampleIDs[str_detect(sampleIDs,"rbf")]
kallistoDirsRBF=kallistoDirs[str_detect(kallistoDirs,"rbf")]

sampleIDsRNA=sampleIDs[str_detect(sampleIDs,"trna")]
kallistoDirsRNA=kallistoDirs[str_detect(kallistoDirs,"trna")]

# 2. handle metadata reading and path association
s2cRBF=read.table(metadataFileRBF,header=TRUE,stringsAsFactors=FALSE,sep=",")
s2cRBF=dplyr::mutate(s2cRBF,path=kallistoDirsRBF)

s2cRNA=read.table(metadataFileRNA,header=TRUE,stringsAsFactors=FALSE,sep=",")
s2cRNA=dplyr::mutate(s2cRNA,path=kallistoDirsRNA)

# 3. running Sleuth
soRBF=sleuth_prep(s2cRBF,extra_bootstrap_summary=TRUE,num_cores=1)

soRNA=sleuth_prep(s2cRNA,extra_bootstrap_summary=TRUE,num_cores=1)
