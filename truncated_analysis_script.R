library(SummarizedExperiment)
library("maxstat")
library("survival")
library(readxl)
library(TCGAbiolinks)
library("writexl")
library(tidyverse)
library(survminer)
library(dplyr)
library(janitor)





################start here? bing in rdataframe from prior analyis
df1 = rownames_to_column(df1, var = "submitter_id")
sample_ids <- read_excel("C:/TCGA/iDEP/protein_coding/protein_coding_genes2.xlsx")

my_second_df <- filter(df1, submitter_id %in% sample_ids$gene_id)

genedata2 = my_second_df

genedata2 = t(genedata2)
genedata3 = as.data.frame(genedata2[-c(2:3),]) # bit at end removes random rows

#This makes the first column a real column 
data4 = rownames_to_column(genedata3, var = "submitter_id")

data4 = janitor::row_to_names(data4, 1, remove_rows_above = FALSE)

####removes unnecessary patient identified text
patient<-gsub("normalized_count_","",as.character(data4$submitter_id))
patient = as.data.frame(patient)
patient = substr(patient$patient,1,nchar(patient$patient)-16)
patient = as.data.frame(patient)

bind = cbind(data4,patient)
bind2 = as.data.frame(bind)


######joins clinical and gene expression data
dff <- cp4 %>% full_join(bind2)

dff_clean = dff[-c(1:6,8) ]

dff_final = as.data.frame(dff_clean)

dff_final$splitMGE = gsub('0','low', dff_final$splitMGE)
dff_final$splitMGE = gsub('1','high', dff_final$splitMGE)


dfftr = t(dff_final)

dfftr = as.data.frame(dfftr)


write.csv(dfftr, "c:/TCGA/iDEP/truncated_files/foridep.csv", row.names = TRUE)



