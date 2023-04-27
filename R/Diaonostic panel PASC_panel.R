#------------------------------------------------
#Load Library
#Loading required packages and setting global options
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(lmerTest))
suppressPackageStartupMessages(library(rstatix))
suppressPackageStartupMessages(library(parallel))

options(stringsAsFactors = FALSE, useFancyQuotes = FALSE)
#------------------------------------------------

############################################################################
#
# read data
# target protein list
#
plist <- read_excel("DiagnosticPanelSelection.xlsx", sheet = 5) %>%
  as.data.frame() 
plist

plist <- plist[1:36,]
names(plist)[1] <-"Protein"
plist

table(plist$commonwithJH)


#
# AIFI metadata
#
adata <- pivotDF %>% as_tibble()
adata

any(duplicated(adata$PTID))
with(adata, table(Group, iGroup, useNA = "ifany"))

adata <- adata %>% mutate(Sample = PTID_visit)
names(adata)
adata
dim(adata)


#
# AIFI Olink data
#

read_txt_file <- function(filename){
  xdata <- fread(filename) %>% t() %>% as.data.frame()
  names(xdata) <- xdata[1,]
  xdata <- xdata[-1,] 
  
  samples <- rownames(xdata)
  
  xdata <- xdata %>% lapply(as.double) %>% as.data.frame()
  xdata$Sample <- samples
  
  return(xdata)
}

# read data
aodata <- t(Olink_rulein) %>%
          as.data.frame() %>%
          rownames_to_column(var = "Sample")
any(duplicated(aodata$Sample))

# proteins of interest
aodata <- aodata %>% select(c("Sample", plist$Protein)) 
aodata
dim(aodata)
names(aodata)


#
# merge
#
dim(adata)
dim(aodata)
adata <- merge(adata, aodata, by="Sample")
dim(adata)
adata

######################################################################
# subset full Olink on selected proteins, trasnpose matric (rows as PTID_visit, columns as proteins, merge with metadata)
######################################################################
#
# Jim Heath metadata
#
jdata <- fread("JH_Olinkmat_ovlp_matannot_col.txt") %>% as_tibble()
jdata

jdata <- jdata %>% mutate(Sample = Subject)
names(jdata)

table(jdata$clusterID3, useNA = "ifany")
with(jdata, table(`Healthy or INCOV`, clusterID3, useNA = "ifany"))
with(jdata, table(symptomatic, clusterID3, useNA = "ifany"))


#
# JH Olink data
#

jodata <- read_txt_file("JH_Olinkmat_full.txt")
any(duplicated(jodata$Sample))

# proteins of interest
names(jodata)
names(jodata) <- names(jodata) %>% tstrsplit("_") %>% .[[1]]
jodata <- jodata[,names(jodata) %in% c("Sample", plist$Protein)]
names(jodata)

# duplicate measurements on IL6
vals <- jodata %>% select(c(IL6, IL6.1, IL6.2)) %>% apply(1, median)
jodata <- jodata %>% mutate(IL6 = vals) %>% select(-c(IL6.1, IL6.2)) 
names(jodata)


#
# merge
#
dim(jdata)
dim(jodata)
jdata <- merge(jdata, jodata, by="Sample")
dim(jdata)
names(jdata)


######################################################################
# subset JH Olink on selected proteins, trasnpose matric (rows as PTID_visit, columns as proteins, merge with metadata)
######################################################################

################################################################################################
#
# inflammatory vs. non-inflammatory
#

#
# common proteins only
#
proteins <- plist$Protein %>% intersect(names(adata)) %>% intersect(names(jdata))
#proteins <- proteins[1:9]
proteins <- proteins[1:15]


###########################################################################
## take 15 common proteins and subset on adata, jdata, add Type (0,1) column to adata and jdata and name it tdata and vdata respectively
###########################################################################
#
# training: AIFI PASC only
#
with(adata, table(Group, iGroup, useNA = "ifany")) 

tdata <- adata %>% filter(Group =="Infected PASC") %>% 
  mutate(Type=ifelse(iGroup=="Inflammatory Group", 1, 0))

dim(tdata)
with(tdata, table(Type, iGroup, useNA = "ifany"))
names(tdata)

tdata <- tdata  %>% select(c(names(tdata)[1:18], "Type", proteins)) 
dim(tdata)
names(tdata)


#
# validation data: JH PASC only
#
with(jdata, table(symptomatic, clusterID3, useNA = "ifany"))

vdata <- jdata %>% filter(symptomatic == "INCOV yes") %>%
  mutate(Type = ifelse(clusterID3=="E", 1, 0))

dim(vdata)
with(vdata, table(Type, clusterID3, useNA = "ifany"))
names(vdata)

vdata <- vdata  %>% select(c(names(vdata)[1:10], "Type", proteins)) 
dim(vdata)
names(vdata)


#
# panel selection
#

library(tictoc)

source("Diaonostic panel LogisticRegression.R")
source("Diaonostic panel Performance.R")


# parameters
#set.seed(20220619)
set.seed(20220810)
nfold <- 10
nboots <- 10000
#nboots <- 100


# forward selection
tic()
xvResults <- xvSelection(tdata, "Type", proteins, nboots, nfold, MC_cores=4)
toc()


plot.XV(xvResults)
summary.XV(xvResults)
write.XV(xvResults, "inflammatory_PASC_panels.csv")



#
# Best model
#
bmodel <- xvResults[[14]] %>% .["Median",colnames(.) != "AUC"]
write.csv(t(bmodel), "inflammatory_PASC_coefficients.csv", row.names=F)

# performance on training 
tdata$Score <- addModelScore(tdata, model=bmodel)

pltROC_scores(tdata, "Type", "Score", "", legends="training", cols="red")
boxplot(Score ~ Type, data=tdata)


# performance on validation 
vdata$Score <- addModelScore(vdata, model=bmodel)

pltROC_scores(vdata, "Type", "Score", "", legends="validation", cols="red")
boxplot(Score ~ Type, data=vdata)
hist(vdata$Score)
