

### Analyse de l'expression différentielle par le logiciel DESeq2


library(RColorBrewer)
library(gplots)
library(limma)
library(FactoMineR)
library(dplyr)
library(ggplot2)
library(DESeq2)
source("https://bioconductor.org/biocLite.R")
library("AnnotationDbi")

library(pathview)
library(gage)
library(gageData)
biocLite("gage")








args = commandArgs(trailingOnly=TRUE)

pathoutput=args[1]

if (args[3]=="Mouse"){

	library("org.Mm.eg.db")
	data(kegg.sets.mm)
	data(sigmet.idx.mm)
	idd="mmu"
	kegg.sets = kegg.sets.mm[sigmet.idx.mm]
	org.eg.db<- org.Mm.eg.db
	print("on est chez la souris")
} else {

	library("org.Hs.eg.db")
	data(kegg.sets.hs)
	data(sigmet.idx.hs)
	idd="hsa"
	kegg.sets = kegg.sets.hs[sigmet.idx.hs]
	org.eg.db<- org.Hs.eg.db
        print("on est chez l'humain")
}


allcounts= read.delim(args[2],sep=" ", as.is=T, check.names=F)

annotations=allcounts[,]
counts = allcounts[,3:ncol(allcounts)]
rownames(counts)=annotations$"GeneId"


## Read in sample information

serie=read.delim(paste0(args[1],"serie.txt"), stringsAsFactors=F)
serie=arrange(serie, condition, replicat)

(condition <- serie$condition)



ind = match(serie$SampleName,colnames(counts))
counts = counts[,ind]


## remove zeros: at least 5 counts in at least 3 samples
keep = which(apply(counts, 1, function(x){sum(x>5)>3})) 
counts = counts[keep,]
annotations = annotations[keep,]

## DESeq2 analysis: Standard
(coldata <- data.frame(row.names=colnames(counts), condition))
dds=DESeqDataSetFromMatrix(countData=counts, colData=coldata, design = ~ condition)

dds = DESeq(dds)
normcounts = counts(dds,normalized=T)
colnames(normcounts) = colnames(counts)


ThePath=pathoutput
write.table(normcounts, file=paste(ThePath,"/Normalized_Counts.xls", sep=""),  sep='\t', quote=F, row.names=T)



initialisation<-function(condition1,condition2){
  print(condition1)
  print(condition2)
  res1 = results(dds, contrast=c("condition",condition1,condition2), 
                 cooksCutoff = FALSE, independentFiltering = TRUE) #, lfcThreshold=0.5) ## 1
  res1=data.frame(res1)
  res1 = data.frame(Symbol = annotations$GeneId,
                    res1)
  
  res1$absFC = abs(res1$log2FoldChange)
  
  
  return(res1)
}

#visual graphics
graphics<-function(res1, condition1, condition2){

  plot(log2(res1$baseMean), res1$log2FoldChange, pch=16, xlab="Log2 baseMean", ylab="Log2 FC", main=paste("MA-plot: ",condition1," vs ",condition2))
  points(log2(res1$baseMean)[res1$padj<0.05], 
         res1$log2FoldChange[res1$padj<0.05],
         pch=16, col=2)
  legend("topright", "Adj-P<0.05", pch=16,col=2)

  plot(res1$log2FoldChange, -log10(res1$padj), pch=16, 
       xlab="Log2 FC", ylab="Log10 pvalue", 
       main=paste("volcano-plot: ",condition1," vs ",condition2))
  points(res1$log2FoldChange[res1$padj<0.05], 
         -log10(res1$padj)[res1$padj<0.05],
         pch=16, col=2)
  legend("topleft", "Adj-P<0.05", pch=16,col=2)

}  
  
  
selection_gene<-function(res1, condition1, condition2){
  # order by BH adjusted p-value
  resOrdered<- res1[order(res1$padj),]
  
  # how many differentially expressed genes ? FDR=5%, |fold-change|>2 (up and down)
  # get differentially expressed gene matrix
  sig <- resOrdered[!is.na(resOrdered$padj) &  resOrdered$padj<0.05 ,] #&  abs(resOrdered$log2FoldChange)>=1,]
  
  
  # select genes
  selected <-   rownames(sig)
  
  x<-log2(counts(dds,normalized=TRUE)[rownames(dds) %in% selected,])
  x[x==-Inf] <- 0 
  
  aaa<-rownames(subset(serie, condition==condition1 | condition==condition2))
  
  heatmap.2(x[,as.numeric(aaa)], scale="row",Rowv = TRUE, Colv= FALSE, dendrogram="row", trace="none", margin=c(4,6), cexRow=0.5, cexCol=0.7, keysize=1 )
  return(selected)
}





pathways<- function(condition1,condition2, res1,selected, kegg.sets, org.eg.db, allcounts){


  ## Pathway

  res2=res1[rownames(res1) %in% selected,]

  colnames(res2)[colnames(res2)=="Symbol"] <- "GeneId"
  file_to_annoted=merge(x = res2, y = allcounts[ , c("GeneId", "gene_name")], by = "GeneId", all.x=TRUE)


 

  #res1$entrez = select(org.Mm.eg.db, keys=row.names(res1),  column="ENTREZID",keytype="SYMBOL", multiVals="first")
  gggg<-select(org.eg.db, keys=file_to_annoted$gene_name,  column="ENTREZID",keytype="SYMBOL", multiVals="first")
  gggg<-gggg[match(unique(gggg$SYMBOL),gggg$SYMBOL),]

  file_to_annoted$entrez=gggg$ENTREZID
  file_to_annoted$name = mapIds(org.eg.db, keys=file_to_annoted$gene_name,  column="GENENAME",keytype="SYMBOL",multiVals="first")
  #file_to_annoted$GO = mapIds(org.eg.db, keys=file_to_annoted$gene_name,  column="GO",keytype="SYMBOL",multiVals="list")
  file_to_annoted$UNIPROT=mapIds(org.eg.db, keys=file_to_annoted$gene_name,  column="UNIPROT",keytype="SYMBOL",multiVals="first")
  file_to_annoted$PATH=as.character(mapIds(org.eg.db, keys=file_to_annoted$gene_name,  column="PATH",keytype="SYMBOL",multiVals="list"))


  print(head(file_to_annoted))
  write.table(file_to_annoted, file=paste("ALL_",condition1,"_vs_",condition2,".xls", sep=""),
              sep="\t", row.names=F)
  


  foldchanges = file_to_annoted$log2FoldChange
  names(foldchanges) = file_to_annoted$entrez
  print(foldchanges)
  
  # Get the results
  keggres = gage(foldchanges, gsets=kegg.sets, same.dir=TRUE)
  
  # Look at both up (greater), down (less), and statistics.
  lapply(keggres, head)
  write.table(keggres$greater, file = "keggres.txt",sep = "\t")
   
  # Get the pathways upregulate
  keggrespathways = data.frame(id=rownames(keggres$greater), keggres$greater) %>% 
    tbl_df() %>% 
    filter(row_number()<=5) %>% 
    .$id %>% 
    as.character()
  keggrespathways
  
  # Get the IDs.
  keggresids = substr(keggrespathways, start=1, stop=8)
  #keggresids
  
 
  
  # plot multiple pathways (plots saved to disk and returns a throwaway list object) !!!!!!!!!!!!!!!!!!!!!!!!!!
	tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, kegg.native =T, same.layer=F, species=idd, min.nnodes = 0, out.suffix=paste(condition1,"_vs_",condition2,"_upregulate_n",match(pid,keggresids), sep="")))
  
  
    
  
  
  # Get the pathways downregulate
  keggrespathways = data.frame(id=rownames(keggres$less), keggres$less) %>% 
    tbl_df() %>% 
    filter(row_number()<=5) %>% 
    .$id %>% 
    as.character()
  keggrespathways
  
  # Get the IDs.
  keggresids = substr(keggrespathways, start=1, stop=8)
  keggresids
  
  
  # plot multiple pathways (plots saved to disk and returns a throwaway list object) !!!!
  tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, kegg.native =T, same.layer=F, species="mmu", min.nnodes = 0, out.suffix=paste(condition1,"_vs_",condition2,"_downregulate_n",match(pid,keggresids), sep="")))

}






  
differentialexpression<-function(condition1,condition2){
  dir.create(paste0(args[1],condition1,"_vs_",condition2), showWarnings = FALSE)
  setwd(paste0(args[1],condition1,"_vs_",condition2))
  pdf(file=paste0(condition1,"_vs_",condition2,"_FoldchangeRplots.pdf"))
  res1=initialisation(condition1,condition2)
  graphics(res1, condition1,condition2)
  selected=selection_gene(res1, condition1, condition2)
  

  
  dat=res1[rownames(res1) %in% selected,]
  
  write.table(res1[rownames(res1) %in% selected,], file=paste(condition1,"_vs_",condition2,".xls", sep=""),
              sep="\t", quote=F, row.names=F)

  dev.off()
  pathways(condition1,condition2, res1, selected, kegg.sets, org.eg.db, allcounts)
  setwd(args[1])
}

conditionToTest=read.delim(paste0(args[1],"compare.txt"), stringsAsFactors=F, sep=" ")


for (i in 1:nrow(conditionToTest)){
differentialexpression(toString(conditionToTest$condition1[i]), toString(conditionToTest$condition2[i]))
}


