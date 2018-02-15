#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(seqTools)


path=getwd()


filenames <- dir(path=args[1], pattern="*_R2.fastq.gz")
#filenames <- file.path(args[1],  c("WT_kidney_2_R1.fastq.gz","WT_kidney_3_R1.fastq.gz","kidney_a11_3_R1.fastq.gz"))
setwd(args[1])
#print(filenames)
len=length(filenames)
fq <- fastqq(filenames, k=9 , probeLabel=filenames)

fileNames(fq)

setwd(path)
setwd(args[2])
pdf(file="Kmerplots.pdf")
#Returns integer vector  of  length nFiles  For  each FASTQ file, the absolute number of containes ’N’ nucleotide entries is given.
#nNnucs(fq)
#
#nFiles(fq)
#nReads(fq)
#maxSeqLen(fq)
#collectTime(fq)
#collectDur(fq)
#slc<-seqLenCount(fq)
#nf<-nucFreq(fq,1)
#nf[1:4,1:len]
#seqLen(fq)
#
#kc<-kmerCount(fq)
#kc[1:len,]

#plotKmerCount(fq,c(1:len), file="kmerplot")
#
#ph<-phred(fq,1)
#ph[25:35,1:15]

#pq<-phredQuantiles(fq,c(0.25,0.5,0.75),1)
#plotNucFreq(fq,len)
# Nucleotide count
#plotNucCount(fq,3) 
# GC content
#gcContent(fq,1)
#
#fqq<-fq[1]

aa<- cbDistMatrix(fq)


hc <- hclust(as.dist(aa))
hcd <- as.dendrogram(hc, lty=2, lwd=2)
op <- par(mar = c(20, 1, 1, 5))

plot(hcd, horiz=FALSE, las=1, edgePar=list(lwd=2, lty=1, col="black", main="k=9 R1"))
dev.off()
