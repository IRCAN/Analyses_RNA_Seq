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
#sfdirs <- file.path(data, c("2_S5/2_S5_transcripts_quant","15_S3/15_S3_transcripts_quant","9_S12/9_S12_transcripts_quant","10_S13/10_S13_transcripts_quant","12_S15/12_S15_transcripts_quant","14_S2/14_S2_transcripts_quant","16_S4/16_S4_transcripts_quant","4_S7/4_S7_transcripts_quant","6_S9/6_S9_transcripts_quant","8_S11/8_S11_transcripts_quant",
#"11_S14/11_S14_transcripts_quant","13_S16/13_S16_transcripts_quant","1_S1/1_S1_transcripts_quant","3_S6/3_S6_transcripts_quant","5_S8/5_S8_transcripts_quant","7_S10/7_S10_transcripts_quant"))

prepare_fish_for_sleuth(sfdirs)

library("sleuth")

base_dir <- pathoutput

sample_id <- dir(file.path(base_dir,""))

s2c <- read.table(file.path(base_dir, "info.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::select(s2c, sample =sample, condition, path)

setwd(pathoutput)
#sf_dirs <- dir()  

#s2c <- dplyr::mutate(s2c, path = sf_dirs)

model <- "~ condition"   	 
# load data and fit the model	 
print(s2c)	 
so <- sleuth_prep(s2c, as.formula(model)) %>%
  sleuth_fit()
models(so)


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

so <- sleuth_wt(so, 'conditionTert10j')
so <- sleuth_wt(so, 'conditionControl10j')
