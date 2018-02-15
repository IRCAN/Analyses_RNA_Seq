

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
library(biomaRt)
library(factoextra)
require("biomaRt")
args = commandArgs(trailingOnly=TRUE)

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y, use="pairwise.complete.obs"))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}

mypanel = function(x,y,...)
{
  points(x,y,pch=16,...)
  abline(0,1,col=2)
}

smooth.pairs = function(...) smoothScatter(..., nrpoints = 1000, add = TRUE)
#Recuperation de l'organisme étudié (souris ou humain), pour définir les paramètres adéquats

#------------------------------------------------------------------------------------------------------------------------------------------#
#si nouvel organisme, mettre references pour kegg ci-dessous, selon la meme structure
#-----souris
if (args[2]=="mm10" | args[2]=="mm9"){
	library("org.Mm.eg.db")
	data(kegg.sets.mm)
	data(sigmet.idx.mm)
	idd="mmu"
	kegg.sets = kegg.sets.mm[sigmet.idx.mm]
	org.eg.db<- org.Mm.eg.db
	specie_dataset="mmusculus_gene_ensembl"
	print("on est chez la souris")
#---zebrafish
} else if (args[2]=="zv10"){
	library("org.Dr.eg.db")

	
	#print(class(kegg.sets))
	#data(kegg.sets.ko)
	#data(sigmet.idx.ko)
	data(kegg.gs)
	kg.dre=kegg.gsets("dre")
	idd="dre"
	kegg.sets = kg.dre$kg.sets[kg.dre$sigmet.idx]
	#print((kegg.sets))
	
	#kegg.sets =kegg.gs
	org.eg.db<- org.Dr.eg.db
	
	specie_dataset="drerio_gene_ensembl"
###
#	data(go.sets.dre)
#	data(go.subs.dre)
#	names(go.subs.dre)
#	go.mf=go.sets.dre[go.subs.dre$MF]
#
###


	print("on est chez le poisson zebre")
#----humain
} else if (args[2]=="hg37" | args[2]=="hg38"){

	library("org.Hs.eg.db")
	data(kegg.sets.hs)
	data(sigmet.idx.hs)
	idd="hsa"
	kegg.sets = kegg.sets.hs[sigmet.idx.hs]
	org.eg.db<- org.Hs.eg.db
	specie_dataset="hsapiens_gene_ensembl"
        print("on est chez l'humain")
}







#-----------------------------------------------------------------------------------------------------------------------------------#
#                         Preparation Analayse Expression differentielle a partir des resultats obtenus par rsem
#-----------------------------------------------------------------------------------------------------------------------------------#

library("tximport")
library("readr")
library("tximportData")	

if (args[3]=="rsem"){
	pathoutput=paste0(args[1],"/DiffExpressRSEM/")

} else if (args[3]=="salmon"){
	pathoutput=paste0(args[1],"/DiffExpressKmer/")
}

s2c <- read.table(file.path(pathoutput, "serie.txt"), header = TRUE, stringsAsFactors=FALSE)
path=getwd()


if (args[3]=="rsem"){
	s2c$path= paste0("/",args[1],"/Genes/",s2c$sample,".isoforms.results")

} else if (args[3]=="salmon"){
	s2c$path= paste0("/",args[1],"/",s2c$sample,"/",s2c$sample,"_transcripts_quant/quant.sf")
}


print(s2c$path)
s3c <- s2c$path
names(s3c)<- s2c$sample


t2g=args[4]
tx2gene=read.table(t2g, header=TRUE)
tx2gene <- dplyr::rename(tx2gene,target_id= Transcrit_Id)# , ens_gene = GeneId) 

if (args[3]=="rsem"){
	txi <- tximport(s3c, type="rsem", tx2gene=tx2gene, txOut=F) 
} else if (args[3]=="salmon"){
	txi <- tximport(s3c, type="salmon", tx2gene=tx2gene, txOut=F)
}
	print(names(txi))
	print(head(txi$counts))
	print(head(txi$length))
	txi$length[txi$length == 0] <- 1
	serie=read.delim(paste0(pathoutput,"serie.txt"), stringsAsFactors=F)
	serie=arrange(serie, condition)
	(condition <- serie$condition)


calculDDS<-function(){
	dds <- DESeqDataSetFromTximport(txi,
		                           colData = s2c,
		                           design = ~condition)

	keep <- which(apply(counts(dds), 1, function(x){sum(x>5)>3}))
	dds <- dds[keep,]

	dds = DESeq(dds)

	return(dds)
	
}

#-----------------------------------------------------------------------------------------------------------------------------------#
#                      utilisation de deseq2, normalisation des donnees
#-----------------------------------------------------------------------------------------------------------------------------------#
normalisation<-function(condition1,condition2, dds){

	
	normcount = counts(dds,normalized=T)
	normcounts=as.data.frame(normcount)

	

	ThePath=pathoutput
	print(head((normcounts)))
		  #keep = which(apply(normcounts, 1, function(x){sum(x>5)>3}))
		  #normcounts = normcounts[keep,]
		  #normcounts$target_id=rownames(normcounts)
		#normcounts$"target_id"=rownames(normcounts)
		references=unique((tx2gene[,1:3]))
		 a=((tstrsplit(rownames(normcounts), "\\.")))
	  	normcounts$"GeneId"=a[[1]]
		print(head(references))
		print(head(normcounts))
		 normcounts2=merge(y = normcounts, x = references[ , c("GeneId","ext_gene")], by = "GeneId", all.y=T)
		 print("normcount ok")

	write.table(unique(normcounts2), file=paste(ThePath,"/Normalized_Counts_genes.xls", sep=""),  sep='\t', quote=F, row.names=F)
	return(normcounts)
}
#-----------------------------------------------------------------------------------------------------------------------------#
#
#------------------------------------------------------------------------------------------------------------------------------#
initialisation<-function(condition1,condition2, dds){
	  print(condition1)
	  print(condition2)
	  res1 = results(dds, contrast=c("condition",condition1,condition2), 
	                 cooksCutoff = FALSE, independentFiltering = FALSE, betaPrior=FALSE) #, lfcThreshold=0.5) ## 1
	
	  res1=data.frame(res1)
	  res1$absFC = abs(res1$log2FoldChange)
	print("ok c'est bon")
  	return(res1)
}


#-----------------------------------------------------------------------------------------------------------------------------#
# graphics: comparaison similarite echantillon meme condition, volcano-plot, MA-plot
#------------------------------------------------------------------------------------------------------------------------------#
#visual graphics
graphics<-function(res1, condition1, condition2, normcounts){
		print("graphiques")
		countst=t((normcounts[,1:ncol(normcounts)-1]))
		counts2=(normcounts[,1:ncol(normcounts)-1])

	###ACP###
	res.pca <- prcomp(countst, scale = TRUE)
	groups <- as.factor(serie$condition)
	#p<-fviz_pca_ind(res.pca)
	
	p<-fviz_pca_ind(res.pca,
            col.ind = groups, # colorer par groupes
            #palette = c("#00AFBB",  "#FC4E07"),
            addEllipses = TRUE, # Ellipse de concentration
            ellipse.type = "confidence",
            legend.title = "Groups",
            repel = TRUE
             )
	plot(p)
	##END ACP###
	
	
	


	pairs(log2(counts2[,serie$condition==condition1]+1), upper.panel=panel.cor,
	      main=condition1)


	## Pairwise comparison: counts,  samples
	pairs(log2(counts2[,serie$condition==condition2]+1), upper.panel=panel.cor,
	      main=condition2)


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
  
#-----------------------------------------------------------------------------------------------------------------------------#
# recuperation des genes differentiallement exprimes
#------------------------------------------------------------------------------------------------------------------------------#
selection_gene<-function(res1, condition1, condition2,dds){
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



#-----------------------------------------------------------------------------------------------------------------------------#
# enrichissements, annotations, kegg pathways
#------------------------------------------------------------------------------------------------------------------------------#

pathways<- function(condition1,condition2, res1,selected, kegg.sets, org.eg.db){
	  

	  ## Pathway

	  res2=res1[rownames(res1) %in% selected,]
	  print(head(res2))

		a=((tstrsplit(rownames(res2), "\\.")))
	  	res2$"GeneId"=a[[1]]
		 #res2$"GeneId"=rownames(res2)
		references=unique((tx2gene[,1:3]))
		print(head(references))
         	file_to_annoted=merge(y = res2, x = references[ , c("GeneId","ext_gene")], by = "GeneId", all.y=T)
		print(head(file_to_annoted))

	  file_to_annoted$entrez<-as.character(mapIds(org.eg.db, keys=as.character(file_to_annoted$ext_gene),  column="ENTREZID",keytype="SYMBOL", multiVals="first"))
	  file_to_annoted$name = as.character(mapIds(org.eg.db, keys=as.character(file_to_annoted$ext_gene),  column="GENENAME",keytype="SYMBOL",multiVals="first"))
	  #file_to_annoted$GO = as.character(mapIds(org.eg.db, keys=file_to_annoted$ext_gene,  column="GO",keytype="SYMBOL",multiVals="list"))
	  file_to_annoted$GOALL = as.character(mapIds(org.eg.db, keys=as.character(file_to_annoted$ext_gene),  column="GOALL",keytype="SYMBOL",multiVals="list"))
	  file_to_annoted$GOALL=gsub("\n","",file_to_annoted$GOALL)

	  file_to_annoted$UNIPROT=as.character(mapIds(org.eg.db, keys=as.character(file_to_annoted$ext_gene),  column="UNIPROT",keytype="SYMBOL",multiVals="first"))
	  file_to_annoted$PATH=as.character(mapIds(org.eg.db, keys=as.character(file_to_annoted$ext_gene),  column="PATH",keytype="SYMBOL",multiVals="list"))
	 

	  write.table(unique(file_to_annoted), file=paste("ALL_genes_",condition1,"_vs_",condition2,".xls", sep=""),sep="\t", row.names=F)
	  
	  foldchanges = file_to_annoted$log2FoldChange
	  names(foldchanges) = file_to_annoted$entrez
	  print(head(foldchanges))
	  ##############

	  # Get the results
	 print("ouuuuuuuuuuuuuuuuuuuuuuuuuuuui")
	  keggres = gage(foldchanges, gsets=kegg.sets, same.dir=TRUE)
	  # Look at both up (greater), down (less), and statistics.
	  lapply(keggres, head)
	  write.table(keggres$less, file = "genes_keggresless.txt",sep = "\t")
	  write.table(keggres$greater, file = "genes_keggresgreater.txt",sep = "\t")
	  write.table(keggres, file = "genes_keggres.txt",sep = "\t")
	  a=tbl_df(keggres$greater)
	  b=which(a$q.val < 0.1) ###change qvalue for kegg pathway
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
			tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, kegg.native =T, same.layer=F, species=idd, min.nnodes = 0, out.suffix=paste("genes_V2",condition1,"_vs_",condition2,"_upregulate_n",match(pid,keggresids), sep="")))
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
			  tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, kegg.native =T, same.layer=F, species=idd, min.nnodes = 0, out.suffix=paste("genes_",condition1,"_vs_",condition2,"_downregulate_n",match(pid,keggresids), sep="")))
	 }
}






  
differentialexpression<-function(condition1,condition2){
	  dir.create(paste0(pathoutput,condition1,"_vs_",condition2), showWarnings = FALSE)
	  setwd(paste0(pathoutput,condition1,"_vs_",condition2))
	  pdf(file=paste0("genes_",condition1,"_vs_",condition2,"_FoldchangeRplots.pdf"))
	  dds=calculDDS()
	  normcounts=normalisation(condition1,condition2, dds)
	  res1=initialisation(condition1,condition2, dds)
	  print("ok c'est bon 2")
	  graphics(res1, condition1,condition2, normcounts)
	  selected=selection_gene(res1, condition1, condition2, dds)
	

	  
		a=((tstrsplit(rownames(res1), "\\.")))
		res1$"GeneId"=a[[1]]
		#res1$target_id=rownames(res1)
		colnames(res1)[which(names(res1) == "GeneId")] <- "GeneId"
		references=unique((tx2gene[,1:3]))
         	res3=unique(merge(y = res1, x = references[ , c("GeneId", "ext_gene")], by = "GeneId", all.y=TRUE))
	

	  write.table(res3, file=paste("genes_No_cutoff_",condition1,"_vs_",condition2,".xls", sep=""),
		      sep="\t", quote=F, row.names=F)

	  write.table(res3[rownames(res1) %in% selected,], file=paste("genes_",condition1,"_vs_",condition2,".xls", sep=""),
		      sep="\t", quote=F, row.names=F)

	  ############################ajout orthologues###############################""
	 # human = useMart("ensembl",dataset="drerio_gene_ensembl")
         # orth = getBM(attributes = c("hsapiens_homolog_canonical_transcript_protein","hsapiens_homolog_subtype","hsapiens_homolog_orthology_type","hsapiens_homolog_perc_id"),
	 # filters="ensembl_gene_id",
	 # values= c("ENSG00000213281","ENSG00000171862"), mart = human)


	###############################################################################################################################
	

	  res2=res1[rownames(res1) %in% selected,]
	  print(head(res2))
	  dev.off()
	  pathways(condition1,condition2, res1, selected, kegg.sets, org.eg.db)
	  setwd(firstpath)
}

firstpath=getwd()
conditionToTest=read.delim(paste0(pathoutput,"compare.txt"), stringsAsFactors=F, sep="	")


for (i in 1:nrow(conditionToTest)){
	differentialexpression(toString(conditionToTest$condition1[i]), toString(conditionToTest$condition2[i]))
}


