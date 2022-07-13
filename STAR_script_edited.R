#These are the libraries that need installing (below). You will need to install them from the packages tab

#tcgabiolinks need code below directly into the console to install 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")


library(SummarizedExperiment)
library("maxstat")
library("survival")
library(readxl)
library(TCGAbiolinks)
library("writexl")
library(tidyverse)
library(survminer)
library(dplyr)

#check all the libraries above are activated

#################################################
#This section is code is to query TCGA tumours (and TARGET or BEATAML)

# Making a TCGA query
query <- GDCquery(project = "TCGA-CHOL",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "STAR - Counts")


# Making a TARGET query
query <- GDCquery(project = "TARGET-ALL-P3",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "STAR - Counts")


# Making a BEATAML1.0-COHORT query
query <- GDCquery(project = "CGCI-BLGSP",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "STAR - Counts")

##########################################
####Technical workarounds for TARGET (ignore this block for TCGA)
#####duplicate workaround for TARGET-nbl tumours with duplicates
query1=query
tmp=query$results[[1]]
tmp=tmp[which(!duplicated(tmp$cases)),]
query1$results[[1]]=tmp

GDCdownload(query1, method = "api", files.per.chunk = 3, 
            directory = "C:/TCGA/star")

dataRNA <- GDCprepare(query = query1, directory = "C:/TCGA/star",
                      save = TRUE, save.filename = "dataRNA.RData")

##########################################
#This downloads your query
##########################################
###Download query putting 3 files per chunk into C:/TCGA  (may need to repeat if download is incomplete)
GDCdownload(query, method = "api", files.per.chunk = 3, 
            directory = "C:/TCGA/star")

#alternative download method if api is not working
GDCdownload(query, method = "client", files.per.chunk = 3, 
            directory = "C:/TCGA/")

#if still not working make sure workspace is in C drive folder rather than onedrive

##########################################


####Preparing a large ranged summarized exp
dataRNA <- GDCprepare(query = query, directory = "C:/TCGA/star",
                      save = TRUE, save.filename = "dataRNA.RData")



###to get LargeSimpleList & make dataframe
df = assays(dataRNA)


#########################if there are memory limits
memory.size() ### Checking your memory size
memory.limit() ## Checking the set limit
memory.limit(size=56000) ### expanding your memory 56000 


#Produces a big DF (all genes vs all cases)
df1 = as.data.frame (df)

#########################
#optional writing to excel (including gene column)
data1 = df1
data1$Gene = row.names(data1)
write_xlsx(data1, "c:/TCGA/mygenedata.xlsx")
########################


####This code finds your gene name(s) identifier
GenelistDF = as.data.frame (rowData(dataRNA))
###can filter DF typing in gene of interest
summarygenelist = GenelistDF %>% filter(gene_name == "DMD")


#selecting specific genes from df1, transposing and making new DF

#This selects only tpm_unstranded data which is in group 4 of df1 (there is optional FPKM and FPKM-UQ)
data_frame_mod <- df1[df1$group==4,]


#this then finds DMD using the ENSG identifier
genedata = data_frame_mod[c("ENSG00000198947.16"),] 
genedata2 = t(genedata) # transposes data
genedata3 = as.data.frame(genedata2[-c(1:2),]) # removes two first rows (not needed)



###This extracts the clinical data for each patient
clinicaldata = as.data.frame(colData(dataRNA))

###########This combines clinical and gene expression data using cbind function
cmb1 <- cbind(clinicaldata, genedata3)

##This selects just primary tumours (removing non-tumour and metastatic samples)
cmb1=cmb1 %>% filter(grepl('Primary', sample_type))

#Optional output to excel
write_xlsx(cmb1, "c:/TCGA/clinicalplusgene.xlsx")


################### This collects survival and gene expression data
MGE = cmb1$`genedata2[-c(1:2), ]`   #gene expression
timedeath = cmb1$days_to_death      #time to death
timelast = cmb1$days_to_last_follow_up  # time to last follow up
cens = cmb1$vital_status            #status (dead or alive)
patient = cmb1$patient              #patient



#for binding columns with survival and gene expression data
forcutpoint = cbind(patient, MGE, timedeath, timelast, cens)
cutpointdf = as.data.frame(forcutpoint)

#optional
write_xlsx(cutpointdf, "c:/TCGA/cutpointdf.xlsx")

library(dplyr)

#putting censored and non-censored data in one column
surv = as.data.frame(coalesce(cutpointdf$timedeath,cutpointdf$timelast))
cutpointdf3 = cbind(surv,cutpointdf)

#renamed fitst column to 'time'
names(cutpointdf3)[names(cutpointdf3)=="coalesce(cutpointdf$timedeath, cutpointdf$timelast)"] = "time"

#converting characters to numerics
cutpointdf3$time <- as.numeric(as.character(cutpointdf3$time))  # Convert time variable to numeric
cutpointdf3$MGE <- as.numeric(as.character(cutpointdf3$MGE))  # Convert MGE to numeric
cutpointdf3$timedeath <- as.numeric(as.character(cutpointdf3$timedeath))  # Convert one variable to numeric
cutpointdf3$timelast <- as.numeric(as.character(cutpointdf3$timelast))  # Convert one variable to numeric

#Recodes dead/alive to 0/1 
cutpointdf3$cens<-recode(cutpointdf3$cens, 'Alive'=0, 'Dead'=1)
################################

library("maxstat")
library("survival")

###Tanya = (na.omit(cutpointdf3)) # removes an rows with NA

#This splits data into high and low expressing groups using survival data. There are different optional methods for 
dichot <- maxstat.test(Surv(time, cens) ~ MGE,
                       data=cutpointdf3, smethod="LogRank",
                       pmethod="condMC", B = 9999)

print(dichot) # this prints the cutpoint value, p-value etc

plot(dichot, xlab="mean gene expression") # This plots a cutpoint graph 




###################################################
### code chunk number 2:produce module for next bit
###################################################
mod <- maxstat.test(Surv(time, cens) ~ MGE,
                    data=cutpointdf3, smethod="LogRank",
                    pmethod="condMC", B = 9999)




###################################################
### code chunk number 3: Survival curves
###################################################
splitMGE <- rep(1, nrow(cutpointdf3))
cp4 <- cbind(cutpointdf3, splitMGE)
cp4$splitMGE[cp4$MGE <= mod$estimate] <- 0


###pretty plotting method
ggsurvplot(
  fit = survfit(Surv(time, cens) ~ splitMGE, data = cp4), 
  xlab = "Survival time in days", 
  ylab = "Overall survival probability")

###########################
#super fancy plotting method (adjust xlim value as needed)

ggsurvplot(
  fit = survfit(Surv(time, cens) ~ splitMGE, data = cp4),             # data used to fit survival curves.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimates of survival curves.
  xlim = c(0,7000),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Survival time in days",   # customize X axis label.
  break.time.by = 500,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
)

###########################
#Statistics section


#median survival time
survfit(Surv(time, cens) ~ splitMGE, data = cp4)

#Comparing survival times between groups 
survdiff(Surv(time, cens) ~ splitMGE, data = cp4)

#calculated p=value (log-Rank)
sd <- survdiff(Surv(time, cens) ~ splitMGE, data = cp4)
1 - pchisq(sd$chisq, length(sd$n) - 1)


#univariate cox regression
coxph(Surv(time, cens) ~ splitMGE, data = cp4)


#hazard ratio
coxph(Surv(time, cens) ~ splitMGE, data = cp4) %>% 
  gtsummary::tbl_regression(exp = TRUE)

################################










