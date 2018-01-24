###
### this script defines differentially expressed transcripts using DESeq2
###

library('DESeq2')
library('stringr')

# 0. user defined variables
tag='trna'

setwd("~/github/aukera/extra/expression/pipeline.star.htseq-count.deseq2")
countsDir="/Volumes/omics4tb/alomana/projects/TLR/data/counts"

DESeqResultsFile=paste("/Volumes/omics4tb/alomana/projects/TLR/data/DESeq2/significance",tag,"csv",sep='.')
DESeqNormalizedCountsFile=paste('/Volumes/omics4tb/alomana/projects/TLR/data/DESeq2/normalizedCounts',tag,'csv',sep='.')

# 1. handle data reading
sampleIDs=dir(file.path(countsDir))

# 1.1. selecting data files and data names
sampleSubset=sampleIDs[str_detect(sampleIDs,tag)]
t1=sampleSubset[str_detect(sampleSubset,"tp.1")]
t4=sampleSubset[str_detect(sampleSubset,"tp.4")]
samples=c(t1,t4)
timeCondition=sapply(strsplit(samples,split='.',fixed=TRUE),function(x) (paste('tp',x[5],sep='.')))

# 2. formatting variables
sampleTable=data.frame(sampleName=samples,fileName=samples,condition=timeCondition)
dds=DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,directory=countsDir,design= ~ condition)

# 3. filtering low abundance transcripts
keep=rowSums(counts(dds)) >= 10
dds=dds[keep,]

# 4. performing differential expression
dds=DESeq(dds)
res=results(dds)
summary(res)

resultsNames(dds)
resLFC=lfcShrink(dds,coef="condition_tp.4_vs_tp.1")

res05=results(dds,alpha=0.05)
summary(res05)

# 5. plotting
plotMA(res)
plotMA(resLFC)

plotCounts(dds, gene=which.min(res$padj),intgroup="condition")

# 6. saving CSV files
resOrdered=res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered),file=DESeqResultsFile)

# 7. exporting normalized counts
values=counts(dds,normalized=TRUE)
write.csv(as.data.frame(values),file=DESeqNormalizedCountsFile)
