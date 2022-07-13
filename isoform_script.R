


# Making a TCGA query
query <- GDCquery(project = "TCGA-CHOL",
                  data.category = "Gene expression",
                  legacy = TRUE,
                  data.type = "Isoform expression quantification")

GDCdownload(query, method = "api", files.per.chunk = 3, 
            directory = "C:/TCGA/isoforms")


#####duplicate workaround for TARGET-nbl tumours with duplicates
query1=query
tmp=query$results[[1]]
tmp=tmp[which(!duplicated(tmp$cases)),]
query1$results[[1]]=tmp


dataRNA <- GDCprepare(query = query1, directory = "C:/TCGA/isoforms",
                      save = TRUE, save.filename = "dataRNA.RData")
#############

#just getting normalised isoform vlaues
New.testSet <- dataRNA[,!grepl("scaled", colnames(dataRNA))]
New.testSet2 <- New.testSet [,!grepl("raw", colnames(New.testSet))]

df1 = New.testSet2

#selecting specific genes from df1, transposing and making new DF
genedata = df1[c("uc004dco.1"),] #dp71b
genedata2 = t(genedata)
genedata3 = as.data.frame(genedata2)

######nor working
genedata4 = str_remove_all(genedata3, "[normalized_count_]")
genedata4 = as.data.frame(genedata4)
#######

#####sort genedata3 out


###Getting clinical data 
clinicalquery <- GDCquery(project = "TCGA-CHOL",
                  data.category = "Clinical",
                  legacy = TRUE,
                  data.type = "Clinical data")

GDCdownload(query, method = "api", files.per.chunk = 3, 
            directory = "C:/TCGA/isoforms")

#####duplicate workaround for TARGET-nbl tumours with duplicates???????????
clinquery1=clinicalquery
tmp=clinicalquery$results[[1]]
tmp=tmp[which(!duplicated(tmp$cases)),]
clinquery1$results[[1]]=tmp


clindataRNA <- GDCprepare(query = clinquery1, directory = "C:/TCGA/isoforms",
                      save = TRUE, save.filename = "dataRNA.RData")


#############




df1 = clindataRNA$clinical_patient_chol

df2 = as.data.frame(df1[-c(1:2),]) # bit at end removes random rows

cmb1 <- cbind(genedata3, df2)


##########

MGE = cmb1$`genedata2[-c(1:2), ]`
timedeath = cmb1$days_to_death
timelast = cmb1$days_to_last_follow_up
cens = cmb1$vital_status
patient = cmb1$patient


###########







