#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(seqTools)
print(args)
setwd(args[2])
filenames <- dir(path=args[1], pattern="*R1.fastq.gz")



fq <- fastqq(file.path(args[1], filenames[-9]), probeLabel=c(filenames[-9]), k=9 )

fileNames(fq)
getK(fq)

#Returns integer vector  of  length nFiles  For  each FASTQ file, the absolute number of containes ’N’ nucleotide entries is given.
nNnucs(fq)
#
nFiles(fq)
nReads(fq)
maxSeqLen(fq)
collectTime(fq)
collectDur(fq)
slc<-seqLenCount(fq)
nf<-nucFreq(fq,1)
nf[1:4,1:16]
seqLen(fq)
probeLabel(fq)= c("10jT1","10jT2","5jT1","5jC1","10jC1","10jC2","5jC2","10jC3","5jC3","5jT2","5jC4","10jT3","10jC4","10jT4","5jT3","5jT4")
#probeLabel(fq)<-1:nFiles(fq)
#
kc<-kmerCount(fq)
kc[1:16,]
plotKmerCount(fq,c(1:16))
#
ph<-phred(fq,1)
ph[25:35,1:15]
pq<-phredQuantiles(fq,c(0.25,0.5,0.75),1)
plotNucFreq(fq,16)
# Nucleotide count
plotNucCount(fq,3) 
# GC content
gcContent(fq,1)
#
fqq<-fq[1]


aa<- cbDistMatrix(fq)


hc <- hclust(as.dist(aa))
hcd <- as.dendrogram(hc, lty=2, lwd=2)
op <- par(mar = c(3, 1, 1, 5))
plot(hcd, horiz=FALSE, las=1, edgePar=list(lwd=2, lty=1, col="black", main="k=9 R2"))
