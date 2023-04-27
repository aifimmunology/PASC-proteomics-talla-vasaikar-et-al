# use Modules to perform ssGSEA and identify PASC clusters

# cleanup & set directors
rm(list=ls())

################################################################################
# Loading required packages and setting global options
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(fgsea))
suppressPackageStartupMessages(library(cowplot))
library("GSEABase")
library("GSVA")
library("ComplexHeatmap")
library("circlize")
################################################################################
#Load Olink data
load("Olink_NPX.Rda")
mat <- mat2021

################################################################################
#Perform ssGSEA using modules
geneset <- getGmt("PathwayModules-54.gmt")
temp <- data.frame(do.call(rbind, strsplit(names(geneset), split = "_")), stringsAsFactors = F)

length(geneset)
mat_exp <- mat
set.seed(20220122)
ssgsea_res <- gsva(data.matrix(mat_exp), geneset, method="ssgsea", abs.ranking=FALSE, min.sz=2, max.sz=2000, ssgsea.norm=T)
boxplot(ssgsea_res)
mat <- ssgsea_res %>% as.data.frame()

#Define samples
ann_sub1 <- ann %>% filter(Group == "Uninfected") %>% as.data.frame()
ann_sub2 <- ann %>% filter(Group == "Infected Recovered" & DaysPSO >= 60) %>%
  group_by(PTID) %>% top_n(1, DaysPSO) %>% as.data.frame() #Last timepoint
ann_sub3 <- ann %>% filter(Group == "Infected PASC" & DaysPSO >= 60) %>%
  group_by(PTID) %>% top_n(-1, DaysPSO) %>% as.data.frame() #First timepoint

# selected samples so far
ann_sub <- rbind(ann_sub1, ann_sub2, ann_sub3)

# summary
table(ann_sub$Group)
table(table(ann_sub$PTID))

table(ann$Group)

# clean up
rm(ann_sub1, ann_sub2, ann_sub3, mat_exp)
################################################################################

################################################################################
#combine metadata with Olink data & standardize 
# proteins
proteins <- row.names(mat)

# samples
samples <- names(mat)

# transpose Olink data & assign PTID_Visit
odata <- mat %>% t() %>% as_tibble() %>%  mutate(PTID_visit = samples)

# commbine ann with odata
mdata <- merge(ann, odata, by="PTID_visit", all=T)
dim(mdata)
table(mdata$Group)

# data of selected samples
sdata <- mdata %>% filter(PTID_visit %in% ann_sub$PTID_visit)
dim(sdata)
table(sdata$Group)

# set Type for donor group
sdata <- sdata %>% mutate(Type=ifelse(grepl("PASC", sdata$Group), 1, 0))
with(sdata, table(Group, Type))

# clean up
rm(odata)
################################################################################

################################################################################
# identify rule-in marker for PASC vs. recovered + uninfected
# identify recovered vs PASC rule-in biomarkers using partial AUC

library(pROC)
library(parallel)
library(tictoc)

# parameters
set.seed(20220122)
conf_level <- 0.9
pauc_range <- c(0.9, 1)

# function for evaluation based on rule-in pAUC
pauc_result <- function(types, values, pauc_range, conf_level){
  pauc0 <- 0.5*(pauc_range[2]-pauc_range[1])^2 # expected value from random classifier
  pauc_ci <- as.numeric(ci.auc(types, values, partial.auc=pauc_range, conf.level=conf_level, boot.n=1000, quiet=T))
  return(data.frame(pAUC=pauc_ci[2], pAUC_low=pauc_ci[1], pAUC_hi=pauc_ci[3], Biomarker=ifelse(pauc_ci[1]>pauc0, "Yes", "No")))
}

# collect performance of all proteins
tic()
marker_results <- mclapply(proteins, function(x){pauc_result(sdata$Type, sdata[,x], pauc_range=pauc_range, conf_level=conf_level)},
                           mc.cores=4) %>% do.call(rbind,.) %>%  cbind(data.frame(Marker=proteins), .) 
toc()

table(marker_results$Biomarker)
biomarkers=marker_results[marker_results$Biomarker %in% "Yes",]$Marker
biomarkers

mcols <- setdiff(colnames(sdata), proteins)
bdata <- sdata[,c(mcols, biomarkers)]
dim(bdata)

tdata <- bdata[,biomarkers] %>% scale() %>% as.data.frame()
row.names(tdata) <- bdata$PTID
tdata <- na.omit(tdata)
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

#Symptoms data
library("readxl")
library("dtplyr")
symptom_data <- read_excel("Symptomology.xlsx", sheet="Sheet #2 - Combo Symptom Groups") %>% lazy_dt() %>% as_tibble()
symptom_data <- data.frame(symptom_data)
row.names(symptom_data) <- symptom_data$PTID_visit
symptom_data <- symptom_data[,12:20]
sdata_FTP <- symptom_data[ann_sub_tdata$PTID_visit,]

ha_col <- HeatmapAnnotation(Age=anno_barplot(ann_sub_tdata$Age), 
                            daysPSO=anno_barplot(ann_sub_tdata$DaysPSO), 
                            df=data.frame(Status=ann_sub_tdata$Group,
                                          Sex=ann_sub_tdata$Sex,
                                          age_bin=ann_sub_tdata$age_bin),
                            col = list(Status = c("Uninfected"="grey","Infected Recovered"="darkgreen","Infected PASC"="magenta"),
                                       Sex = c("Male"="skyblue", "Female"="pink")))
ha_col1 <- HeatmapAnnotation(Age=anno_barplot(ann_sub_tdata$Age),
                             daysPSO=anno_barplot(ann_sub_tdata$DaysPSO), 
                             df=data.frame(Status=ann_sub_tdata$Group,
                                           Sex=ann_sub_tdata$Sex,
                                           age_bin=ann_sub_tdata$age_bin,
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
                                        age_bin = c("Young"="orange", "Old"="brown"),
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

set.seed(20220127)               
htmap1 <- Heatmap(t(tdata), cluster_rows = T,  cluster_columns = T,
               row_km = 4, row_km_repeats = 10, column_km = 5, column_km_repeats = 10,
               # clustering_method_columns = "average",
               # column_split = factor(ann_sub$Group, levels=c("Uninfected","Infected Recovered","Infected PASC")),
               na_col = "grey", col = colorRamp2(c(-2,0,2), c("purple","black","yellow")),
               row_names_max_width=unit(10, "cm"), #rect_gp = gpar(col = "white"),
               show_column_names = T, 
               column_names_gp = gpar(fontsize = 5), row_names_gp = gpar(fontsize = 6),
               top_annotation = ha_col1,
               heatmap_legend_param = list(title = "NPX(Z-score)",heatmap_legend_side = "right") )

htmap1

set.seed(20220127)               
rO <- row_order(htmap1)
rO <- colnames(tdata[,as.numeric(unlist(rO))])
set.seed(20220127)                
cO <- column_order(htmap1)
cO <- row.names(tdata[as.numeric(unlist(cO)),])
################################################################################

################################################################################
#Longitudinal data
tdata_mat <- mat[biomarkers,] %>% t() %>% scale() %>% as.data.frame()
ov <- intersect(row.names(tdata_mat), row.names(ann))
ann_LTP <- ann[ov,]
tdata_mat <- tdata_mat[ov,]
symptom_tdata <- symptom_data[ov,]

set.seed(20220127)               
ha_col2 <- HeatmapAnnotation(Age=anno_barplot(ann_LTP$Age), 
                             daysPSO=anno_barplot(ann_LTP$DaysPSO), 
                             df=data.frame(PTID=ann_LTP$PTID,
                                           Status=ann_LTP$Group, Sex=ann_LTP$Sex,
                                           age_bin=ann_LTP$age_bin,
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
                                        age_bin = c("Young"="orange", "Old"="brown"),
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


#by unsupervised
set.seed(20220127) 
htmap2 <- Heatmap(t(tdata_mat), cluster_rows = T,  cluster_columns = T,
                  row_km = 4, row_km_repeats = 10, column_km = 5, column_km_repeats = 10,
                  #clustering_method_columns = "average",
                  #column_split = factor(ann_LTP$PTID, levels = cO),
                  na_col = "grey", col = colorRamp2(c(-2,0,2), c("purple","black","yellow")),
                  row_names_max_width=unit(10, "cm"), #rect_gp = gpar(col = "white"),
                  show_column_names = T, 
                  column_names_gp = gpar(fontsize = 3), row_names_gp = gpar(fontsize = 6),
                  top_annotation = ha_col2,
                  heatmap_legend_param = list(title = "NPX(Z-score)",heatmap_legend_side = "right") )

#by donor
htmap3 <- Heatmap(t(tdata_mat), cluster_rows = F,  cluster_columns = F,
                  #row_km = 4, row_km_repeats = 10, column_km = 5, column_km_repeats = 10,
                  #clustering_method_columns = "average",
                  column_split = factor(ann_LTP$PTID, levels = cO),
                  na_col = "grey", col = colorRamp2(c(-2,0,2), c("purple","black","yellow")),
                  row_names_max_width=unit(10, "cm"), #rect_gp = gpar(col = "white"),
                  show_column_names = T, 
                  column_names_gp = gpar(fontsize = 4), row_names_gp = gpar(fontsize = 6),
                  top_annotation = ha_col2,
                  heatmap_legend_param = list(title = "NPX(Z-score)",heatmap_legend_side = "right") )


pdf("temp3.pdf", width=15, height=10)
print(htmap1)
print(htmap2)
print(htmap3)
dev.off()
################################################################################

################################################################################
#FTP and longitudinal clustering
#Get the cluster information
c_data <- read_excel("c2-Clustering-FTP-Longitudinal.xlsx", sheet="Sheet1") %>% lazy_dt() %>% as_tibble()
L_data <- read_excel("c2-Clustering-FTP-Longitudinal.xlsx", sheet="Sheet2") %>% lazy_dt() %>% as_tibble()
row.names(c_data) <- c_data$PTID_visit
row.names(L_data) <- L_data$PTID_visit
tdata_mat <- mat[biomarkers,] %>% t() %>% scale() %>% as.data.frame()
ov <- intersect(row.names(tdata_mat), row.names(ann))
ann_LTP <- ann[ov,]
tdata_mat <- tdata_mat[ov,]
symptom_tdata <- symptom_data[ov,]
c_data <- c_data[ov,]
L_data <- L_data[ov,]

set.seed(20220127)               
ha_col2 <- HeatmapAnnotation(Age=anno_barplot(ann_LTP$Age), 
                             daysPSO=anno_barplot(ann_LTP$DaysPSO), 
                             df=data.frame(PTID=ann_LTP$PTID,
                                           Status=ann_LTP$Group, Sex=ann_LTP$Sex,
                                           age_bin=ann_LTP$age_bin,
                                           FTP_Cluster=c_data$Cluster1,
                                           Long_Cluster=L_data$Cluster2,
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
                                        age_bin = c("Young"="orange", "Old"="brown"),
                                        FTP_Cluster=c("C1"="darkgreen","C2"="green","C3"="pink","C4"="red","C5"="brown"),
                                        Long_Cluster=c("C1"="darkgreen","C2"="green","C3"="pink","C4"="red","C5"="brown"),
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
#by donor
htmap4 <- Heatmap(t(tdata_mat), cluster_rows = F,  cluster_columns = F,
                  #row_km = 4, row_km_repeats = 10, column_km = 5, column_km_repeats = 10,
                  #clustering_method_columns = "average",
                  column_split = factor(ann_LTP$PTID, levels = cO),
                  na_col = "grey", col = colorRamp2(c(-2,0,2), c("purple","black","yellow")),
                  row_names_max_width=unit(10, "cm"), #rect_gp = gpar(col = "white"),
                  show_column_names = T, 
                  column_names_gp = gpar(fontsize = 4), row_names_gp = gpar(fontsize = 6),
                  top_annotation = ha_col2,
                  heatmap_legend_param = list(title = "NPX(Z-score)",heatmap_legend_side = "right") )
pdf("temp3c.pdf", width=30, height=10)
print(htmap4)
dev.off()
################################################################################

################################################################################
library("readxl")
library("dtplyr")
c_data <- read_excel("c2-Clustering-FTP-Longitudinal.xlsx", sheet="Sheet1") %>% lazy_dt() %>% as_tibble()
L_data <- read_excel("c2-Clustering-FTP-Longitudinal.xlsx", sheet="Sheet2") %>% lazy_dt() %>% as_tibble()

c_data <- c_data %>% data.frame()
L_data <- L_data %>% data.frame()
row.names(c_data) <- c_data$PTID_visit
row.names(L_data) <- L_data$PTID_visit

#temp <- ann_sub_tdata[as.character(c_data$Sample),]
#temp <- ann_LTP[as.character(L_data$Sample),]
#write.csv(temp, file="temp.csv")

table(c_data$Cluster1, c_data$Group)
count <- with(c_data, table(Cluster1, Group))
freq <- t(t(count)/colSums(count))
colSums(freq)
count
freq

#count <- with(co_cluster, table(Cluster, as.character(NSymptoms)))
count <- with(c_data, table(Cluster1, Pulmonary))
count <- with(c_data, table(Cluster1, Fatigue...malaise))
freq <- t(t(count)/colSums(count))
colSums(freq)
count
freq

df <- merge(c_data, L_data, by="PTID_visit")
row.names(df) <- df$PTID.x
df <- df[as.character(c_data$Sample),]
count <- with(df, table(Cluster1, Cluster2))
count
freq <- t(t(count)/colSums(count))
colSums(freq)
freq

#Longitudinal data
#tdata_mat <- mat[biomarkers,] %>% t() %>% scale() %>% as.data.frame()
all.equal(row.names(tdata_mat), row.names(ann_LTP))
c_data_sub <- c_data[row.names(tdata_mat),]
L_data_sub <- L_data[row.names(tdata_mat),]
nModules <- colnames(tdata_mat)

#Module data
outNet <- read_excel("c2-Clustering-FTP-Longitudinal.xlsx", sheet="Sheet3") %>% lazy_dt() %>% as_tibble()
outNet <- outNet %>% data.frame()
row.names(outNet) <- as.character(outNet$module)
outNet <- outNet[nModules,]

splots <- list()
for(i in 1:length(nModules)) {
  nM <- nModules[i]
  temp <- data.frame(sample=row.names(tdata_mat), exp=tdata_mat[,nM],Cluster=c_data_sub$Cluster1, Cluster2=L_data_sub$Cluster2, ann_LTP, stringsAsFactors = F)
  temp1 <- temp[!is.na(temp$Cluster),]
  p1 <- ggplot(temp1, aes(x=Cluster, y=exp)) + geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color=Group), width=0.2) +
    labs(x="", y="ssGSEA score", title=outNet[nM,]$name) +
    theme_bw()
  splots[[i]] <- p1
}

pdf("temp4.pdf", width=30, height=30)
print(plot_grid(plotlist = splots, ncol=5))
dev.off()
################################################################################
