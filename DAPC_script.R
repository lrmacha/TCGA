BRCAquery <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - FPKM-UQ")



####to get large ranged summarized exp
BRCAdataRNA <- GDCprepare(query = BRCAquery, directory = "C:/TCGA",
                      save = TRUE, save.filename = "BRCAdataRNA.RData")
##########################################################################
KIRPquery <- GDCquery(project = "TCGA-KIRP",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - FPKM-UQ")



####to get large ranged summarized exp
KIRPdataRNA <- GDCprepare(query = KIRPquery, directory = "C:/TCGA",
                          save = TRUE, save.filename = "KIRPdataRNA.RData")
##########################################################################

LAMLquery <- GDCquery(project = "TCGA-LAML",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - FPKM-UQ")



####to get large ranged summarized exp
LAMLdataRNA <- GDCprepare(query = LAMLquery, directory = "C:/TCGA",
                          save = TRUE, save.filename = "LAMLdataRNA.RData")
###########################################################################
LGGquery <- GDCquery(project = "TCGA-LGG",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - FPKM-UQ")



####to get large ranged summarized exp
LGGdataRNA <- GDCprepare(query = LGGquery, directory = "C:/TCGA",
                          save = TRUE, save.filename = "LGGdataRNA.RData")
###########################################################################
LUADquery <- GDCquery(project = "TCGA-LUAD",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - FPKM-UQ")



####to get large ranged summarized exp
LUADdataRNA <- GDCprepare(query = LUADquery, directory = "C:/TCGA",
                          save = TRUE, save.filename = "LUADdataRNA.RData")
###########################################################################
PAADquery <- GDCquery(project = "TCGA-PAAD",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - FPKM-UQ")



####to get large ranged summarized exp
PAADdataRNA <- GDCprepare(query = PAADquery, directory = "C:/TCGA",
                          save = TRUE, save.filename = "PAADdataRNA.RData")
###########################################################################
READquery <- GDCquery(project = "TCGA-READ",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - FPKM-UQ")



####to get large ranged summarized exp
READdataRNA <- GDCprepare(query = READquery, directory = "C:/TCGA",
                          save = TRUE, save.filename = "READdataRNA.RData")
############################################################################
THYMquery <- GDCquery(project = "TCGA-THYM",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - FPKM-UQ")



####to get large ranged summarized exp
THYMdataRNA <- GDCprepare(query = THYMquery, directory = "C:/TCGA",
                          save = TRUE, save.filename = "THYMdataRNA.RData")
##########################################################################

UVMquery <- GDCquery(project = "TCGA-UVM",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - FPKM-UQ")



####to get large ranged summarized exp
UVMdataRNA <- GDCprepare(query = UVMquery, directory = "C:/TCGA",
                          save = TRUE, save.filename = "UVMdataRNA.RData")
##########################################################################

