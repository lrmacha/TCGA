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

#########################
#optional writing to excel (including gene column)
data1 = df1
data1$Gene = row.names(data1)
write_xlsx(data1, "c:/TCGA/mygenedata.xlsx")
########################


####to find gene name(s)
GenelistDF = as.data.frame (rowData(dataRNA))
###can filter DF typing in gene of interest
summarygenelist = GenelistDF %>% filter(external_gene_name == "DMD")
print(summarygenelist)


#selecting specific genes from df1, transposing and making new DF
genedata = df1[c("ENSG00000198947"),] 
genedata2 = t(genedata)
genedata3 = as.data.frame(genedata2[-c(1:2),]) # bit at end removes random rows



###extracting clinical data
clinicaldata = as.data.frame(colData(dataRNA))

###########Combining clinical and gene expression data
cmb1 <- cbind(clinicaldata, genedata3)

##removing non-tumour samples
cmb1=cmb1 %>% filter(grepl('Primary', sample_type))

#Optional output to excel
write_xlsx(cmb1, "c:/TCGA/clinicalplusgene.xlsx")
###################
MGE = cmb1$`genedata2[-c(1:2), ]`
timedeath = cmb1$days_to_death
timelast = cmb1$days_to_last_follow_up
cens = cmb1$vital_status
patient = cmb1$patient

#binding columns
forcutpoint = cbind(patient, MGE, timedeath, timelast, cens)
cutpointdf = as.data.frame(forcutpoint)
#optional
write_xlsx(cutpointdf, "c:/TCGA/cutpointdf.xlsx")

library(dplyr)

surv = as.data.frame(coalesce(cutpointdf$timedeath,cutpointdf$timelast))
cutpointdf3 = cbind(surv,cutpointdf)

#renamed column to time
names(cutpointdf3)[names(cutpointdf3)=="coalesce(cutpointdf$timedeath, cutpointdf$timelast)"] = "time"

#converting characters to numerics
cutpointdf3$time <- as.numeric(as.character(cutpointdf3$time))  # Convert time variable to numeric
cutpointdf3$MGE <- as.numeric(as.character(cutpointdf3$MGE))  # Convert MGE to numeric
cutpointdf3$timedeath <- as.numeric(as.character(cutpointdf3$timedeath))  # Convert one variable to numeric
cutpointdf3$timelast <- as.numeric(as.character(cutpointdf3$timelast))  # Convert one variable to numeric

#Recodes dead/alive to 0/1 
cutpointdf3$cens<-recode(cutpointdf3$cens, 'Alive'=0, 'Dead'=1)


library("maxstat")
library("survival")

###Tanya = (na.omit(cutpointdf3)) # removes an rows with NA

dichot <- maxstat.test(Surv(time, cens) ~ MGE,
                         data=cutpointdf3, smethod="LogRank",
                         pmethod="condMC", B = 9999)

print(dichot)

plot(dichot, xlab="mean gene expression")



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

plot(survfit(Surv(time, cens) ~ splitMGE, data=cp4), 
     xlab = "Survival time in days",
     ylab="Probability", cex.lab=1.3, cex.axis=1.3, lwd=2, mark.time = TRUE)
text(80, 0.9, expression("Low DMD gene expression" > 4524), cex=1.3)   
text(80, 0.1, expression("High DMD gene expression" <= 4524 ), cex=1.3)

#https://www.emilyzabor.com/tutorials/survival_analysis_in_r_tutorial.html

survdiff(Surv(time, cens) ~ splitMGE, cp4)

###pretty plotting method
ggsurvplot(
  fit = survfit(Surv(time, cens) ~ splitMGE, data = cp4), 
  xlab = "Survival time in days", 
  ylab = "Overall survival probability")

###########################
ggsurvplot(
  fit = survfit(Surv(time, cens) ~ splitMGE, data = cp4),             # data used to fit survival curves.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimates of survival curves.
  xlim = c(0,6000),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Survival time in days",   # customize X axis label.
  break.time.by = 500,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
)


####super fancy plot
###########################

#median survival time
survfit(Surv(time, cens) ~ splitMGE, data = cp4)

#Comparing survival times between groups 
survdiff(Surv(time, cens) ~ splitMGE, data = cp4)

#calculated p=value
sd <- survdiff(Surv(time, cens) ~ splitMGE, data = cp4)
1 - pchisq(sd$chisq, length(sd$n) - 1)


#cox regression
coxph(Surv(time, cens) ~ splitMGE, data = cp4)

#install summary package

#hazard ratio
coxph(Surv(time, cens) ~ splitMGE, data = cp4) %>% 
  gtsummary::tbl_regression(exp = TRUE)












