library(SummarizedExperiment)
library("maxstat")
library("survival")
library(readxl)
library(TCGAbiolinks)
library("writexl")
library(tidyverse)
library(survminer)
library(dplyr)


# Making a query
query <- GDCquery(project = "TCGA-CHOL",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - FPKM-UQ")


###Download query (may need to repeat if download is incomplete)
GDCdownload(query, method = "api", files.per.chunk = 3, 
            directory = "C:/TCGA/")


####to get large ranged summarized exp
dataRNA <- GDCprepare(query = query, directory = "C:/TCGA",
                      save = TRUE, save.filename = "dataRNA.RData")



###to get LargeSimpleList & make dataframe
df = assays(dataRNA)


#########################if there are memory limits
memory.size() ### Checking your memory size
memory.limit() ## Checking the set limit
memory.limit(size=56000) ### expanding your memory 56000 


#produces DF (all genes vs all cases)
df1 = as.data.frame (df)



################start here? bing in rdataframe from prior analyis
genedata2 = t(df1)
genedata3 = as.data.frame(genedata2[-c(1:2),]) # bit at end removes random rows

#This makes the first column a real column 
data4 = rownames_to_column(genedata3, var = "submitter_id")

####removes unnecessary patient identified text
patient<-gsub("normalized_count_","",as.character(data4$submitter_id))
patient = as.data.frame(patient)
patient = substr(patient$patient,1,nchar(patient$patient)-16)
patient = as.data.frame(patient)


bind = cbind(data4,patient)
bind2 = as.data.frame(bind)



######joins clinical and gene expression data
dff <- cp4 %>% full_join(bind2)

dfftr = t(dff)

dff_clean = dfftr[-c(1:6,8), ]

dff_final = as.data.frame(dff_clean)

#This makes the first column a real column 
dff_final = rownames_to_column(dff_final, var = " ")


write_xlsx(dff_final, "c:/TCGA/foridep.xlsx")

write.table(dff_final, file = "c:/TCGA/foridep.txt", sep = "")




