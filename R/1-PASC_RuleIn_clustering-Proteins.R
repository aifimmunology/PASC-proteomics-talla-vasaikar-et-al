# RuleIn criteria used to select proteins as biomarkers. These biomarkers used
# to cluster PASC data for cluster separation

# cleanup & set directors
rm(list=ls())

################################################################################
# Loading required packages and setting global options
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(ggplot2))

options(stringsAsFactors = FALSE, useFancyQuotes = FALSE, warn = -1)
################################################################################


# load olink data
load("data/Olink_NPX.Rda")
mat <- mat2021

################################################################################
# identify rule-in marker for PASC vs. recovered + uninfected
# identify recovered vs PASC rule-in biomarkers using partial AUC

library(pROC)
library(parallel)
library(tictoc)

# parameters
set.seed(20220122)
conf_level <- 0.99
pauc_range <- c(0.9, 1)

# function for evaluation based on rule-in pAUC
pauc_result <- function(types, values, pauc_range, conf_level){
  pauc0 <- 0.5*(pauc_range[2]-pauc_range[1])^2 # expected value from random classifier
  pauc_ci <- as.numeric(ci.auc(types, values, partial.auc=pauc_range, conf.level=conf_level, boot.n=10000, quiet=T))
  return(data.frame(pAUC=pauc_ci[2], pAUC_low=pauc_ci[1], pAUC_hi=pauc_ci[3], Biomarker=ifelse(pauc_ci[1]>pauc0, "Yes", "No")))
}

# collect performance of all proteins
tic()
marker_results <- mclapply(proteins, function(x){pauc_result(sdata$Type, sdata[,x], pauc_range=pauc_range, conf_level=conf_level)},
                           mc.cores=4) %>%
  do.call(rbind,.) %>% 
  cbind(data.frame(Marker=proteins), .) 
toc()

table(marker_results$Biomarker)
#write.csv(marker_results, "PASC_vs_Others_marker_performance.csv", row.names = F)
################################################################################

################################################################################
# unsupervised clustering of all donors based on biomarkers of PASC vs. recovered + uninfected

# biomarker data
#biomarkers <- marker_results %>% filter(Biomarker=="Yes") %>% select(Marker) %>% unlist() %>% as.vector()
marker_results <- read.csv("data/PASC_vs_Others_marker_performance.csv", stringsAsFactors = F)
biomarkers=marker_results[marker_results$Biomarker %in% "Yes",]$Marker
biomarkers

mcols <- setdiff(colnames(sdata), proteins)
bdata <- sdata[,c(mcols, biomarkers)]
dim(bdata)

# unsupervised clustering
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization

# compute gap statistic
set.seed(20220123)
tdata <- bdata[,biomarkers] %>% scale() %>% as.data.frame()
row.names(tdata) <- bdata$PTID

tdata <- na.omit(tdata)
k21 <- kmeans(tdata, centers = 4, nstart = 100)
fviz_cluster(k21, data = tdata)

kdata <- data.frame(PTID=names(k21$cluster), Cluster=k21$cluster)

kdata <- merge(bdata[,mcols], kdata, by="PTID", all.y=T)
count <- with(kdata, table(Cluster, Group))
freq <- t(t(count)/colSums(count))
colSums(freq)
count
freq

fviz_nbclust(tdata, kmeans, method = "wss")
fviz_nbclust(tdata, kmeans, method = "silhouette")

fviz_nbclust(t(tdata), kmeans, method = "wss")
fviz_nbclust(t(tdata), kmeans, method = "silhouette")

gap_stat <- clusGap(tdata, FUN = kmeans, nstart = 100,
                    K.max = 10, B = 50)
fviz_gap_stat(gap_stat)
################################################################################

################################################################################
# heatmap
library("ComplexHeatmap")
library("circlize")

all.equal(as.character(ann_sub$PTID), row.names(tdata))
ann_sub_tdata <- ann_sub
row.names(ann_sub_tdata) <- ann_sub_tdata$PTID
ann_sub_tdata <- ann_sub_tdata[row.names(tdata),]
all.equal(as.character(ann_sub_tdata$PTID), row.names(tdata))

ha_col <- HeatmapAnnotation(df=data.frame(Status=ann_sub_tdata$Group, Sex=ann_sub_tdata$Sex),
                            col = list(Status = c("Uninfected"="grey","Infected Recovered"="darkgreen","Infected PASC"="magenta"),
                                       Sex = c("Male"="skyblue", "Female"="pink")))

set.seed(20220123)                 
htmap1 <- Heatmap(t(tdata), cluster_rows = T,  cluster_columns = T,
                 row_km = 4, row_km_repeats = 10, column_km = 5, column_km_repeats = 10,
                 # clustering_method_columns = "average",
               # column_split = factor(ann_sub$Group, levels=c("Uninfected","Infected Recovered","Infected PASC")),
               na_col = "grey", col = colorRamp2(c(-2,0,2), c("purple","black","yellow")),
               row_names_max_width=unit(10, "cm"), #rect_gp = gpar(col = "white"),
               show_column_names = T, 
               column_names_gp = gpar(fontsize = 5), row_names_gp = gpar(fontsize = 6),
               top_annotation = ha_col,
               heatmap_legend_param = list(title = "NPX(Z-score)",heatmap_legend_side = "right") )

htmap1

set.seed(20220123)                
rO <- row_order(htmap1)
rO <- colnames(tdata[,as.numeric(unlist(rO))])
set.seed(20220123)                
cO <- column_order(htmap1)
cO <- row.names(tdata[as.numeric(unlist(cO)),])

#Symptoms data
library("readxl")
library("dtplyr")
symptom_data <- read_excel("data/Symptomology.xlsx", sheet="Sheet #2 - Combo Symptom Groups") %>% lazy_dt() %>% as_tibble()
symptom_data <- data.frame(symptom_data)
row.names(symptom_data) <- symptom_data$PTID_visit
symptom_data <- symptom_data[,12:20]
sdata_FTP <- symptom_data[ann_sub_tdata$PTID_visit,]

#FTP annotation
ha_col2 <- HeatmapAnnotation(df=data.frame(Status=ann_sub_tdata$Group, Sex=ann_sub_tdata$Sex,
                                           #Fever_chills=sdata_FTP$Fever...Chills,
                                           #ApetiteLoss=sdata_FTP$Loss.of.Appetite.or.Weight.loss,
                                           Fatigue=sdata_FTP$Fatigue...malaise,
                                           Pulmonary=sdata_FTP$Pulmonary,
                                           Cardiovascular=sdata_FTP$Cardiovascular,
                                           Gastrointestinal=sdata_FTP$Gastrointestinal,
                                           Musculoskeletal=sdata_FTP$Musculoskeletal,
                                           Neurologic=sdata_FTP$Neurologic,
                                           Any.mild.symptom=sdata_FTP$Any.mild.symptom),
                             col = list(Status = c("Uninfected"="grey","Infected Recovered"="darkgreen","Infected PASC"="magenta"),
                                        Sex = c("Male"="skyblue", "Female"="pink"),
                                        Fever_chills=c("Yes"="black", "No"="white"),
                                        ApetiteLoss=c("Yes"="black", "No"="white"),
                                        Fatigue=c("Yes"="black", "No"="white"),
                                        Pulmonary=c("Yes"="black", "No"="white"),
                                        Cardiovascular=c("Yes"="black", "No"="white"),
                                        Gastrointestinal=c("Yes"="black", "No"="white"),
                                        Musculoskeletal=c("Yes"="black", "No"="white"),
                                        Neurologic=c("Yes"="black", "No"="white"),
                                        Any.mild.symptom=c("Yes"="black", "No"="white")
                             ),
                             gp = gpar(col = "black"))

set.seed(20220123)                
htmap2 <- Heatmap(t(tdata), cluster_rows = T,  cluster_columns = T,
                 row_km = 4, row_km_repeats = 10, column_km = 4, column_km_repeats = 10,
                 clustering_method_columns = "average",
                 # column_split = factor(ann_sub$Group, levels=c("Uninfected","Infected Recovered","Infected PASC")),
                 na_col = "grey", col = colorRamp2(c(-2,0,2), c("purple","black","yellow")),
                 row_names_max_width=unit(10, "cm"), #rect_gp = gpar(col = "white"),
                 show_column_names = T, 
                 column_names_gp = gpar(fontsize = 5), row_names_gp = gpar(fontsize = 6),
                 top_annotation = ha_col2,
                 heatmap_legend_param = list(title = "NPX(Z-score)",heatmap_legend_side = "right") )

htmap2
################################################################################

################################################################################
#longitudinal clustering
# Get the protein matrix data
tdata_mat <- mat[biomarkers,] %>% t() %>% scale() %>% as.data.frame()
tdata_mat <- tdata_mat[,rO] #same row order

ov <- intersect(row.names(tdata_mat), row.names(ann))
ann_LTP <- ann[ov,]
tdata_mat <- tdata_mat[ov,]
symptom_tdata <- symptom_data[ov,]

ha_col3 <- HeatmapAnnotation(df=data.frame(PTID=ann_LTP$PTID, Status=ann_LTP$Group, Sex=ann_LTP$Sex,
                                           #Fever_chills=symptom_tdata$Fever...Chills,
                                           #ApetiteLoss=symptom_tdata$Loss.of.Appetite.or.Weight.loss,
                                           Fatigue=symptom_tdata$Fatigue...malaise,
                                           Pulmonary=symptom_tdata$Pulmonary,
                                           Cardiovascular=symptom_tdata$Cardiovascular,
                                           Gastrointestinal=symptom_tdata$Gastrointestinal,
                                           Musculoskeletal=symptom_tdata$Musculoskeletal,
                                           Neurologic=symptom_tdata$Neurologic,
                                           Any.mild.symptom=symptom_tdata$Any.mild.symptom),
                             col = list(Status = c("Uninfected"="grey","Infected Recovered"="darkgreen","Infected PASC"="magenta"),
                                        Sex = c("Male"="skyblue", "Female"="pink"),
                                        Fever_chills=c("Yes"="black", "No"="white"),
                                        ApetiteLoss=c("Yes"="black", "No"="white"),
                                        Fatigue=c("Yes"="black", "No"="white"),
                                        Pulmonary=c("Yes"="black", "No"="white"),
                                        Cardiovascular=c("Yes"="black", "No"="white"),
                                        Gastrointestinal=c("Yes"="black", "No"="white"),
                                        Musculoskeletal=c("Yes"="black", "No"="white"),
                                        Neurologic=c("Yes"="black", "No"="white"),
                                        Any.mild.symptom=c("Yes"="black", "No"="white")
                             ),
                             gp = gpar(col = "black"))

htmap3 <- Heatmap(t(tdata_mat), cluster_rows = T,  cluster_columns = T,
                  row_km = 4, row_km_repeats = 10, column_km = 5, column_km_repeats = 10,
                  clustering_method_columns = "average",
                  #column_split = factor(ann_LTP$PTID, levels = cO),
                  na_col = "grey", col = colorRamp2(c(-2,0,2), c("purple","black","yellow")),
                  row_names_max_width=unit(10, "cm"), #rect_gp = gpar(col = "white"),
                  show_column_names = T, 
                  column_names_gp = gpar(fontsize = 5), row_names_gp = gpar(fontsize = 6),
                  top_annotation = ha_col3,
                  heatmap_legend_param = list(title = "NPX(scaled)",heatmap_legend_side = "right") )


pdf("htmap3.pdf", width=30, height=10)
print(htmap3)
dev.off()
################################################################################
