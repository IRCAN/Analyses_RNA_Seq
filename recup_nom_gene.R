	


args = commandArgs(trailingOnly=TRUE)

if (args[3]=="mm10" | args[3]=="mm9"){
	specie_dataset="mmusculus_gene_ensembl"
} else {
	specie_dataset="hsapiens_gene_ensembl"
}


file_to_annote=read.table(args[1], header=TRUE)


mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
		                 dataset = specie_dataset,
		                 host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
		                             "external_gene_name"), mart = mart)

tx2gene <- dplyr::rename(t2g, Transcrit_Id = ensembl_transcript_id,
		             GeneId = ensembl_gene_id, ext_gene = external_gene_name)

file_annoted=try(merge( x = tx2gene[ , c("Transcrit_Id","ext_gene")], y = file_to_annote,by = "Transcrit_Id", all.y=TRUE))
print(class(file_annoted))
if (class(file_annoted) == "try-error") {
    
	file_annoted=merge( x = tx2gene[ , c("GeneId","ext_gene")], y = file_to_annote,by = "GeneId", all.y=TRUE)
}
write.table(unique(file_annoted), file=paste(args[2]), sep="\t", row.names=F, quote=FALSE)
