

### Analyse de l'expression différentielle par le logiciel DESeq2

library(RColorBrewer)
library(gplots)
library(limma)
library(FactoMineR)
library(dplyr)
library(ggplot2)
library(DESeq2)
library("AnnotationDbi")
library(pathview)
library(gage)
library(gageData)
library("data.table")


args = commandArgs(trailingOnly=TRUE)


#Recuperation de l'organisme étudié (souris ou humain), pour définir les paramètres adéquats
if (args[3]=="mm10" | args[3]=="mm9"){
	library("org.Mm.eg.db")
	data(kegg.sets.mm)
	data(sigmet.idx.mm)
	idd="mmu"
	kegg.sets = kegg.sets.mm[sigmet.idx.mm]
	org.eg.db<- org.Mm.eg.db
	specie_dataset="mmusculus_gene_ensembl"
	print("on est chez la souris")
} else {

	library("org.Hs.eg.db")
	data(kegg.sets.hs)
	data(sigmet.idx.hs)
	idd="hsa"
	kegg.sets = kegg.sets.hs[sigmet.idx.hs]
	org.eg.db<- org.Hs.eg.db
	specie_dataset="hsapiens_gene_ensembl"
        print("on est chez l'humain")
}



if (args[4]=="salmon"){
	pathoutput=paste0(args[1],"/DiffExpressSalmon/")
	library("tximport")
	library("readr")
	library("tximportData")

	

	
	s2c <- read.table(file.path(pathoutput, "serie.txt"), header = TRUE, stringsAsFactors=FALSE)
	path=getwd()
	s2c$path= paste0(path,"/",args[1],"/",s2c$sample,"/",s2c$sample,"_transcripts_quant/quant.sf")

	print(s2c$path)
	s3c <- s2c$path
	names(s3c)<- s2c$sample
	
	allcounts= read.delim(args[2],sep="\t", as.is=T, check.names=F)
	print(head(allcounts))
	counts = allcounts[,3:ncol(allcounts)]
	





	mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
		                 dataset = specie_dataset,
		                 host = 'ensembl.org')


	t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name", "transcript_version"), mart = mart)
	t2g <- within(t2g, target_id <- paste(ensembl_transcript_id, transcript_version,sep='.'))
	t2g <- dplyr::rename(t2g, ens_gene = ensembl_gene_id, ext_gene = external_gene_name) 
	tx2gene <- t2g[,c(5,2,3)]
	#write.table(tx2gene, file=paste("jdshdjs.xls", sep=""),
	#	      sep="\t", quote=F, row.names=T)
	txi <- tximport(s3c, type="salmon", tx2gene=tx2gene, txOut=TRUE)
	serie=read.delim(paste0(pathoutput,"serie.txt"), stringsAsFactors=F)
	serie=arrange(serie, condition)

	(condition <- serie$condition)
	
	dds <- DESeqDataSetFromTximport(txi,
		                           colData = s2c,
		                           design = ~condition)



} else {

	pathoutput=paste0(args[1],"/DiffExpress/")
	allcounts= read.delim(args[2],sep="", as.is=T, check.names=F)
	print(head(allcounts[,2:ncol(allcounts)]))
	annotations=allcounts[,]
	counts = allcounts[,2:ncol(allcounts)]
	rownames(counts)=annotations$"GeneId"

	allcounts2= read.delim(args[5],sep="\t", as.is=T, check.names=F)
	references=(allcounts2[,1:2])
	## Read in sample information
	
	serie=read.delim(paste0(pathoutput,"serie.txt"), stringsAsFactors=F)
	serie=arrange(serie, condition)

	(condition <- serie$condition)

	ind = match(serie$sample,colnames(counts))
	counts = counts[,ind]

	## remove zeros: at least 5 counts in at least 3 samples
	keep = which(apply(counts, 1, function(x){sum(x>5)>3})) 
	counts = counts[keep,]
	annotations = annotations[keep,]
	## DESeq2 analysis: Standard
	(coldata <- data.frame(row.names=colnames(counts), condition))
	dds=DESeqDataSetFromMatrix(countData=counts, colData=coldata, design = ~ condition)
}




dds = DESeq(dds)
normcount = counts(dds,normalized=T)
normcounts=as.data.frame(normcount)
colnames(normcounts) = colnames(counts)
ThePath=pathoutput
print(head(rownames(normcounts)))
if (args[4]=="salmon"){
	  normcounts$target_id=rownames(normcounts)
	  normcounts$target_id2=(substr(normcounts$target_id,1,18))
	  tx2gene$target_id2=substr(tx2gene$target_id,1,18)
	  normcounts2=merge(y = normcounts, x = tx2gene[ , c("target_id2","ext_gene")], by = "target_id2", all.y=TRUE)
}else{
	  
	  a=((tstrsplit(rownames(normcounts), "\\.")))
	  normcounts$"GeneId"=a[[1]]
	  normcounts2=merge(y = normcounts, x = references[ , c("GeneId", "ext_gene")], by = "GeneId", all.y=TRUE)
	  }	
write.table(normcounts2, file=paste(ThePath,"/Normalized_Counts.xls", sep=""),  sep='\t', quote=F, row.names=F)



initialisation<-function(condition1,condition2){
	print(condition1)
	  print(condition2)
	  res1 = results(dds, contrast=c("condition",condition1,condition2), 
	                 cooksCutoff = FALSE, independentFiltering = FALSE, betaPrior=FALSE) #, lfcThreshold=0.5) ## 1
	  res1=data.frame(res1)
	  #res1 = data.frame(Symbol = annotations$GeneId,
	  #                  res1)
  
	  res1$absFC = abs(res1$log2FoldChange)
  
  
  	return(res1)
}

#visual graphics
graphics<-function(res1, condition1, condition2){

  	plot(log2(res1$baseMean), res1$log2FoldChange, pch=16, xlab="Log2 baseMean", ylab="Log2 FC", main=paste("MA-plot: ",condition1," vs ",condition2))
	  points(log2(res1$baseMean)[res1$padj<0.05], 
       	 res1$log2FoldChange[res1$padj<0.05],pch=16, col=2)
  	legend("topright", "padj<0.05", pch=16,col=2)

  plot(res1$log2FoldChange, -log10(res1$padj), pch=16, 
       xlab="Log2 FC", ylab="Log10 padj", 
       main=paste("volcano-plot: ",condition1," vs ",condition2))
  	points(res1$log2FoldChange[res1$padj<0.05], 
         -log10(res1$padj)[res1$padj<0.05],
         pch=16, col=2)
  	legend("topleft", "padj<0.05", pch=16,col=2)

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

	  if (length(x)>2){
		try(heatmap.2(x[,as.numeric(aaa)], scale="row",Rowv = TRUE, Colv= FALSE, dendrogram="row", trace="none", margin=c(4,6), cexRow=0.5, cexCol=0.7, keysize=1 ))
	  }
	  return(selected)
}





pathways<- function(condition1,condition2, res1,selected, kegg.sets, org.eg.db, allcounts){


	  ## Pathway
	
	  res2=res1[rownames(res1) %in% selected,]
	  if (args[4]=="salmon"){
		  res2$target_id=rownames(res2)
		  res2$target_id2=(substr(res2$target_id,1,18))
		  tx2gene$target_id2=substr(tx2gene$target_id,1,18)
		  print(head(res2))
		  print(head(tx2gene))
		  file_to_annoted=merge(y = res2, x = tx2gene[ , c("target_id2","ext_gene")], by = "target_id2", all.y=TRUE)
	  }else{

	 	  a=((tstrsplit(rownames(res2), "\\.")))
		  res2$"GeneId"=a[[1]]
		  file_to_annoted=merge(y = res2, x = references[ , c("GeneId", "ext_gene")], by = "GeneId", all.y=TRUE)
	  }	  
	  #res1$entrez = select(org.Mm.eg.db, keys=row.names(res1),  column="ENTREZID",keytype="SYMBOL", multiVals="first")
	  file_to_annoted$entrez<-as.character(mapIds(org.eg.db, keys=file_to_annoted$ext_gene,  column="ENTREZID",keytype="SYMBOL", multiVals="first"))
	  #gggg<-gggg[match(unique(gggg$SYMBOL),gggg$SYMBOL),]
	  #file_to_annoted$entrez=gggg$ENTREZID
	  file_to_annoted$name = as.character(mapIds(org.eg.db, keys=file_to_annoted$ext_gene,  column="GENENAME",keytype="SYMBOL",multiVals="first"))
	  file_to_annoted$GO = as.character(mapIds(org.eg.db, keys=file_to_annoted$ext_gene,  column="GO",keytype="SYMBOL",multiVals="list"))
	  
	  
	  file_to_annoted$UNIPROT=as.character(mapIds(org.eg.db, keys=file_to_annoted$ext_gene,  column="UNIPROT",keytype="SYMBOL",multiVals="first"))
	  file_to_annoted$PATH=as.character(mapIds(org.eg.db, keys=file_to_annoted$ext_gene,  column="PATH",keytype="SYMBOL",multiVals="list"))

	  file_to_annoted$GO=gsub("\n","",file_to_annoted$GO)
	  
	  write.table(file_to_annoted, file=paste("ALL_",condition1,"_vs_",condition2,".xls", sep=""),sep="\t", row.names=F)
	  

	  

	  foldchanges = file_to_annoted$log2FoldChange
	  names(foldchanges) = file_to_annoted$entrez

	  
	  # Get the results
	  keggres = gage(foldchanges, gsets=kegg.sets, same.dir=TRUE)
	  
	  # Look at both up (greater), down (less), and statistics.
	  lapply(keggres, head)
	  write.table(keggres$less, file = "keggresless.txt",sep = "\t")
	  write.table(keggres$greater, file = "keggresgreater.txt",sep = "\t")
	  write.table(keggres, file = "keggres.txt",sep = "\t")
	  a=tbl_df(keggres$greater)
	  b=which(a$q.val < 0.1)
	  keggresId=FALSE
	 if (length(b)>0){
		  keggresId=TRUE
		  # Get the pathways upregulate
		  keggrespathways = data.frame(id=rownames(keggres$greater), keggres$greater) %>% 
		    tbl_df() %>% 
		    filter(row_number()<=length(b)) %>% 
		    .$id %>% 
		    as.character()
		  keggrespathways
		  
		  # Get the IDs.
		  keggresids = substr(keggrespathways, start=1, stop=8)
		  #keggresids
	  }
	 
	  if (keggresId==TRUE){
		  # plot multiple pathways (plots saved to disk and returns a throwaway list object) !!!!!!!!!!!!!!!!!!!!!!!!!!
			tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, kegg.native =T, same.layer=F, species=idd, min.nnodes = 0, out.suffix=paste(condition1,"_vs_",condition2,"_upregulate_n",match(pid,keggresids), sep="")))
	  }
	  keggresId=FALSE
	    
	  a=tbl_df(keggres$less)
	  b=which(a$q.val < 0.1)
	  if (length(b)>0){
		  keggresId=TRUE
		  # Get the pathways downregulate
		  keggrespathways = data.frame(id=rownames(keggres$less), keggres$less) %>% 
		    tbl_df() %>% 
		    filter(row_number()<=length(b)) %>% 
		    .$id %>% 
		    as.character()
		  keggrespathways
		  
		  # Get the IDs.
		  keggresids = substr(keggrespathways, start=1, stop=8)
		  keggresids
	  }
	  if (keggresId==TRUE){
			  # plot multiple pathways (plots saved to disk and returns a throwaway list object) !!!!
			  tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, kegg.native =T, same.layer=F, species="mmu", min.nnodes = 0, out.suffix=paste(condition1,"_vs_",condition2,"_downregulate_n",match(pid,keggresids), sep="")))
	 }
}






  
differentialexpression<-function(condition1,condition2){
	  dir.create(paste0(pathoutput,condition1,"_vs_",condition2), showWarnings = FALSE)
	  setwd(paste0(pathoutput,condition1,"_vs_",condition2))
	  pdf(file=paste0(condition1,"_vs_",condition2,"_FoldchangeRplots.pdf"))
	  res1=initialisation(condition1,condition2)
	  graphics(res1, condition1,condition2)
	  selected=selection_gene(res1, condition1, condition2)


	  
	  dat=res1[rownames(res1) %in% selected,]
	  if (args[4]=="salmon"){
		  res1$target_id=rownames(res1)
		  res1$target_id2=(substr(res1$target_id,1,18))
		  tx2gene$target_id2=substr(tx2gene$target_id,1,18)
		  res3=merge(y = res1, x = tx2gene[ , c("target_id2","ext_gene")], by = "target_id2", all.y=TRUE)
	  }else{
		  a=((tstrsplit(rownames(res1), "\\.")))
		  res1$"GeneId"=a[[1]]

		  res3=merge(y = res1, x = references[ , c("GeneId", "ext_gene")], by = "GeneId", all.y=TRUE)
	  }	  
	  
	  write.table(res3, file=paste("No_selection_",condition1,"_vs_",condition2,".xls", sep=""),
		      sep="\t", quote=F, row.names=F)

	  write.table(res3[rownames(res1) %in% selected,], file=paste(condition1,"_vs_",condition2,".xls", sep=""),
		      sep="\t", quote=F, row.names=F)

	  dev.off()
	  pathways(condition1,condition2, res1, selected, kegg.sets, org.eg.db, allcounts)
	  setwd(firstpath)
}


firstpath=getwd()
conditionToTest=read.delim(paste0(pathoutput,"compare.txt"), stringsAsFactors=F, sep="	")


for (i in 1:nrow(conditionToTest)){
	differentialexpression(toString(conditionToTest$condition1[i]), toString(conditionToTest$condition2[i]))
}


