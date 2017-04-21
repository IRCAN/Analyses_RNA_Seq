
source("https://bioconductor.org/biocLite.R")
library("AnnotationDbi")

library(pathview)
library(gage)
library(gageData)
biocLite("gage")
library(dplyr)





if args[2]=="Mouse":
	library("org.Mm.eg.db")
	data(kegg.sets.mm)
	data(sigmet.idx.mm)
	kegg.sets = kegg.sets.mm[sigmet.idx.mm]
	org.eg.db<- org.Mm.eg.db
else:
	library("org.Hs.eg.db")
	data(kegg.sets.hs)
	data(sigmet.idx.hs)
	kegg.sets = kegg.sets.hs[sigmet.idx.hs]
	org.eg.db<- org.Hs.eg.db

allcounts= read.delim(args[2],sep=" ", as.is=T, check.names=F)
res2=allcounts2


annotations<-function(kegg.sets, res2, org.eg.db, allcounts)

	
	#res1$entrez = select(org.Mm.eg.db, keys=row.names(res1),  column="ENTREZID",keytype="SYMBOL", multiVals="first")
	gggg<-select(org.eg.db, keys=as.character(res2$GeneID),  column="SYMBOL",keytype="SYMBOL", multiVals="first")
	gggg<-gggg[match(unique(gggg$SYMBOL),gggg$SYMBOL),]
	#head(gggg)
	res2$Symbol=gggg$SYMBOL
	res2$name = mapIds(org.eg.db, keys=res2$Symbol,  column="GENENAME",keytype="SYMBOL",multiVals="first")
	res2$GO = mapIds(org.eg.db, keys=res2$Symbol,  column="GO",keytype="SYMBOL",multiVals="first")
	res2$UNIPROT=mapIds(org.eg.db, keys=res2$Symbol,  column="UNIPROT",keytype="SYMBOL",multiVals="first")
	res2$PATH=mapIds(org.eg.db, keys=res2$Symbol,  column="PATH",keytype="SYMBOL",multiVals="list")
	
	foldchanges = allcounts$FoldChange
	names(foldchanges) = allcounts$GeneID

	# Get the results
	keggres = gage(foldchanges, gsets=kegg.sets, same.dir=TRUE)

	# Look at both up (greater), down (less), and statistics.
	lapply(keggres, head)
	write.table(keggres$greater, file = "keggres.txt",sep = "\t")

	write.table(rbind(keggres$greater, keggres$less), file = "keggres2.txt", sep = "\t")
	keggres.kegg.sig<-sigGeneSet(keggres, outname="keggres.kegg")
	write.table(keggres.kegg.sig, file = "keggres3.txt", sep = "\t")

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

	subset(res2,row.names(res2)%in%(grep("04974",res2$PATH))
	       # Define plotting function for applying later
	       #plot_pathway = function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="mmu", new.signature=FALSE)
	       
	       # plot multiple pathways (plots saved to disk and returns a throwaway list object)
	       tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, limit = list(gene=c(-1,1)), both.dir=F, pathway.id=pid, kegg.native =T, same.layer=F, species="mmu", out.suffix=paste("Upregulate_n",match(pid,keggresids), sep="")))
	       tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, limit = list(gene=1), pathway.id=pid, kegg.native =T, same.layer=F, species="mmu", out.suffix=paste("Upregulate_",match(pid,keggresids),"PV_", sep="")))
	       
	       tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, limit = list(gene=3), pathway.id=pid, kegg.native =F, split.group =T, species="mmu", out.suffix=paste("Splitgroup_upregulate_n",match(pid,keggresids), sep="")))
	       
	       tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, limit = list(gene=3), color=c(1:3), pathway.id=pid, kegg.native =F, expand.node =T, species="mmu", out.suffix=paste("Expandnode_upregulate_n",match(pid,keggresids), sep="")))
	       
	       
	       
	       
	       
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
	       
