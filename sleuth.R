library(wasabi)


args = commandArgs(trailingOnly=TRUE)

pathoutput=args[1]
listesample=args[4]
print(listesample)
if (args[3]=="mm10" | args[3]=="mm9"){
	specie_dataset="mmusculus_gene_ensembl"
} else {
specie_dataset="hsapiens_gene_ensembl"
}
a=(as.list(strsplit(listesample, ","))[[1]])
print(a)
data <- pathoutput
sfdirs <- file.path(data, c(a))

prepare_fish_for_sleuth(sfdirs)

library("sleuth")

base_dir <- pathoutput

sample_id <- dir(file.path(base_dir,""))

s2c <- read.table(file.path(base_dir, "/DiffExpressSalmon/serie.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::select(s2c, sample =sample, condition)

path=getwd()
s2c$path= paste0(path,"/",args[1],"/",s2c$sample,"/",s2c$sample,"_transcripts_quant/")
setwd(pathoutput)
#sf_dirs <- dir()  

#s2c <- dplyr::mutate(s2c, path = sf_dirs)

model <- "~ condition"   	 
# load data and fit the model	 
print(s2c)
so <- sleuth_prep(s2c, ~condition)
print("ok") 
so <- sleuth_prep(s2c, as.formula(model)) %>%
  sleuth_fit()
models(so)
 print ("ouioi")

mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = specie_dataset,
                         host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
so <- sleuth_prep(s2c, ~ condition, target_mapping = t2g)

so <- sleuth_fit(so)
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
results_table <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
sleuth_live(so)

so <- sleuth_wt(so, 'control')

