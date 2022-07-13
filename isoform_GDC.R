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
library(DT)
library(tibble)
library(ggplot2)


query <- GDCquery(project = "TCGA-HNSC",
                  legacy = TRUE,
                  data.category = "Gene expression" ,
                  data.type = "Isoform expression quantification",
                  sample.type = c("Primary Tumor") ,
                  file.type  = c("normalized"))

#####Breast only...there were two RNAseq options
query <- GDCquery(project = "TCGA-BRCA",
                  legacy = TRUE,
                  data.category = "Gene expression" ,
                  data.type = "Isoform expression quantification",
                  sample.type = c("Primary Tumor") ,
                  file.type  = c("normalized") ,
                  experimental.strategy = "RNA-Seq")
#######



#####LAML only
query <- GDCquery(project = "TCGA-LAML",
                  legacy = TRUE,
                  data.category = "Gene expression" ,
                  data.type = "Isoform expression quantification",
                  file.type  = c("normalized"),
                  platform = "Illumina HiSeq")
#########
                  

###Download query (may need to repeat if download is incomplete)
GDCdownload(query, method = "api", files.per.chunk = 3, 
            directory = "C:/TCGA/isoforms")

####to get large ranged summarized exp
dataRNA <- GDCprepare(query = query, directory = "C:/TCGA/isoforms",
                      save = TRUE, save.filename = "dataRNA.RData", )



###to get LargeSimpleList & make dataframe
##df = assays(dataRNA)

df1 = as.data.frame (dataRNA)

########################getting clinical data (change cancer!!!!!)

clinical <- GDCquery_clinic(project = "TCGA-HNSC", 
                            type = "clinical")



####gene expresssion of isoforms
genedata = df1[c("uc004dcy.1","uc004dda.1","uc004dcz.2","uc004ddb.1","uc004dcw.2","uc004dcx.2","uc004dct.1","uc004dcu.1","uc004dcv.1","uc004dcr.1","uc004dcs.1","uc004dcq.1","uc004dcm.1","uc004dcn.1","uc004dco.1","uc004dcp.1","uc011mkb.1","uc004ddd.1","uc010ngp.1","uc010ngr.1"
),] 
genedata2 = t(genedata)
genedata3 = as.data.frame(genedata2)
write_xlsx(genedata3, "c:/TCGA/isogeneexpression.xlsx")

#########


#selecting specific genes from df1, transposing and making new DF
genedata = df1[c("uc004dda.1"),] #example uc004dcm.1 = dp71
genedata2 = t(genedata)
genedata3 = as.data.frame(genedata2)


########
# uc004dcy.1	Dp427p1,	mRNA.
# uc004dda.1	Dp427m,	mRNA.
# uc004dcz.2	Dp427p2,	mRNA.
# uc004ddb.1	Dp427c,	mRNA.
# uc004dcw.2	Dp260-2,	mRNA.
# uc004dcx.2	Dp260-1,	mRNA.
# uc004dct.1	Dp140,	mRNA.
# uc004dcu.1	Dp140b,	mRNA.
# uc004dcv.1	D140ab,	mRNA.
# uc004dcr.1	Dp140c,	mRNA.
# uc004dcs.1	Dp140bc,	mRNA.
# uc004dcq.1	Dp116,	mRNA.
# uc004dcm.1	Dp71,	mRNA.
# uc004dcn.1	Dp71a,	mRNA.
# uc004dco.1	Dp71b,	mRNA.
# uc004dcp.1	Dp71ab,	mRNA.
# uc011mkb.1	Dp40,	mRNA.
# uc004ddd.1		
# uc010ngp.1		
# uc010ngr.1		
##########

###########tidying up text

#This makes the first column a real column 
data4 = rownames_to_column(genedata3, var = "submitter_id")

####removes unnecessary patient identified text
play<-gsub("normalized_count_","",as.character(data4$submitter_id))
play = as.data.frame(play)
play = substr(play$play,1,nchar(play$play)-16)
play2 = as.data.frame(play)


bind = cbind(data4,play2)
bind2 = as.data.frame(bind [,2:3])


######for expression analysis only
bind = cbind(genedata3,clinical)
bind2 = as.data.frame(bind [,2:22])
#######



bind3 = bind2 %>% 
  rename(
    submitter_id = play)

######joins clinical and gene expression data
dff <- clinical %>% full_join(bind3)

##could I remame last column to MGE
names(dff)[length(names(dff))]<-"MGE" 


####change isotype
dff2 = dff[!is.na(dff$MGE),] # removed rows with no gene expression data 
cmb1 = dff2



###################################################################

###################
MGE = cmb1$MGE #needs changing for each isoform
timedeath = cmb1$days_to_death
timelast = cmb1$days_to_last_follow_up
cens = cmb1$vital_status
patient = cmb1$submitter_id



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

#https://www.emilyzabor.com/tutorials/survival_analysis_in_r_tutorial.html


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
  xlim = c(0,8000),         # present narrower X axis, but not affect
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


##### write out cp4 file for graphpad
write_xlsx(cp4, "c:/TCGA/cp4.xlsx")




                          