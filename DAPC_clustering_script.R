library(readxl)
DAPC_HR_heatmap <- read_excel("C:/Users/lrmacha/OneDrive - The University of Nottingham/Module 5 Research project/Thesis_Lee/DAPC_HR_heatmap.xlsx", 
                              sheet = "Heatmap")


df = DAPC_HR_heatmap


##############works for clusting tumours based on DAPC marker set
df_trans = t(df)
df_trans = as.data.frame(df_trans)

hc = df_trans [-1,]
hc <- hclust(dist(hc, method="euclidean"), method="ward.D2")


fviz_dend(hc, k = 3, # Cut in three groups
          cex = 1.0, # label size
          k_colors = "jco",
          color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE) # Add rectangle around groups


##other clustering options
#hc, method="euclidean"), method="single"))
#hc, method="euclidean"), method="complete"))
#hc, method="minkowski", p=1/4), method="ward.D2"))

#colors
#k_colors = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A"))

##Allowed values include also "grey" for grey color palettes; brewer
##palettes e.g. "RdBu", "Blues", ...; and scientific journal palettes from ggsci R
##package, e.g.: "npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons"
##and "rickandmorty".

################



##############works for clusting markers based on tumour set
hc = column_to_rownames(df, var = "...1")

hc <- hclust(dist(hc, method="euclidean"), method="ward.D2")


fviz_dend(hc, k = 3, # Cut in three groups
          cex = 0.9, # label size
          #k_colors = NULL,
          #k_colors = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A"),
          k_colors = "jco",
          color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE) # Add rectangle around groups


##other clustering options
#hc, method="euclidean"), method="single"))
#hc, method="euclidean"), method="complete"))
#hc, method="minkowski", p=1/4), method="ward.D2"))

#colors
#k_colors = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A"))

##Allowed values include also "grey" for grey color palettes; brewer
##palettes e.g. "RdBu", "Blues", ...; and scientific journal palettes from ggsci R
##package, e.g.: "npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons"
##and "rickandmorty".


################

