# Loading required packages and setting global options
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(car))
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(quantmod))
suppressPackageStartupMessages(library(package = BiocParallel))
suppressPackageStartupMessages(library(package = fgsea))
suppressPackageStartupMessages(library(package = GSVA))
suppressPackageStartupMessages(library(package = igraph))

options(stringsAsFactors = FALSE, useFancyQuotes = FALSE, warn = -1)

# make pivot sample metadata from of the rule-in timepoint samples
pivotDF <- read_excel("c2-Clustering-FTP-Longitudinal.xlsx", sheet = 1) %>%
           as.data.frame()

# read full metadata table
tableS1 <- read_excel("Table S1.xlsx", sheet = 1) %>% as.data.frame()



## Clustering based on symptomology : only showing PASC
pivotDF2 <- read_excel("c2-Clustering-FTP-Longitudinal.xlsx", sheet = 2) %>%
            as.data.frame() %>%
            filter(Group == "Infected PASC" & daysbin == ">=60 Days PSO")
symptomsDF <- pivotDF2[, c(2,11:20)] %>%
              gather(symptom, value, -PTID) %>%
              group_by(PTID, symptom) %>%
              mutate(value2 = ifelse("Yes" %in% value, "Yes", value)) %>%
              as.data.frame() %>%
              select(PTID, symptom, value2) %>%
              unique() %>%
              spread(symptom, value2)

symptomsDF[symptomsDF == "Yes"] <- 1
symptomsDF[symptomsDF == "No"] <- 0
symptomsDF[symptomsDF == "NA"] <- 0
symptomsDF[, c(2:ncol(symptomsDF))] <- apply(symptomsDF[,
                                                        c(2:ncol(symptomsDF))],
                                                        2,
                                                        as.numeric)
symptomsDF <- symptomsDF %>% column_to_rownames(var = "PTID")
symptomsDF <- as.data.frame(t(symptomsDF))

# column annotation for heatmap
matannot_col <- pivotDF2 %>% select(PTID, Sex) %>% unique()
rownames(matannot_col) <- NULL
matannot_col <- matannot_col %>% column_to_rownames(var = "PTID")
matannot_col <- matannot_col[colnames(symptomsDF), , drop = F]
table(rownames(matannot_col) == colnames(symptomsDF))

# annotation colors
ann_colors <- list(Sex = c("Male" = "lightblue", "Female" = "lightpink"))

# define custom color palatte
colorPalette <- c("white", "black")

# plot
outFile <- "Fig S2A.pdf"
pdf(file = outFile, width = 15, height = 10)
 pheatmap(mat = symptomsDF,
              color = colorPalette,
              cellwidth = 5,
              cellheight = 10,
              cluster_cols = TRUE,
              cluster_rows = TRUE,
              show_rownames = TRUE,
              show_colnames = TRUE,
              annotation_col = matannot_col,
              annotation_colors = ann_colors,
              fontsize_col = 6,
              fontsize_row = 10,
              fontsize_number = 20,
              border_color = "lightgrey")
dev.off()



####### stacked bar plot of number of patients in groups per cluster
nsamplecluster <- as.data.frame(table(pivotDF$Cluster))
plotDF <- as.data.frame(table(pivotDF$Cluster, pivotDF$Group)) %>%
          mutate(N = nsamplecluster$Freq[match(Var1,
                                        table = nsamplecluster$Var1)],
                 prop = round((Freq/N)*100, 2))


# arrange by cluster with most to least proportion of PASC
porder <- plotDF %>%
          filter(Var2 == "Infected PASC") %>%
          arrange(desc(prop)) %>%
          .$Var1

fCols <- c("Uninfected" = "grey",
           "Infected Recovered" = "black",
           "Infected PASC" = "magenta")

plotbar <- ggplot(data = plotDF,
                   mapping =  aes(x = Var1,
                                  y = prop,
                                  fill = Var2)) +
            scale_fill_manual(values = fCols) +
            geom_bar(position = "stack", stat = "identity") +
            scale_x_discrete(limits = rev(porder)) +
            labs(x = "Cluster number",
                 y = "Proportion of participants per group",
                 fill = "Group") +
            theme_bw() +
            theme(axis.text.x = element_text(size = 12, color = "black"),
                  axis.text.y = element_text(size = 12, color = "black"),
                  axis.line = element_line(color = "grey"),
                  panel.background = element_blank(),
                  panel.grid = element_blank())
pdf(file = "Figure S3A.pdf", width = 5, height = 4)
print(plotbar)
dev.off()


### ICS
ics <- read_excel("ICS/ICS.xlsx", sheet = 1) %>%
       as.data.frame() %>%
       filter(PTID_visit %in% pivotDF$PTID_visit) %>%
       filter(!stimulation == "NA" & Group == "Infected PASC") %>%
       mutate(clusterID = pivotDF$Cluster[match(PTID_visit,
                                    table = pivotDF$PTID_visit)],
              pctpos_adj_CD4 = as.numeric(pctpos_adj_CD4),
              pctpos_adj_CD8 = as.numeric(pctpos_adj_CD8))
ics$stimulation <- factor(ics$stimulation,
                            levels = c("Spike", "Non-Spike",
                                       "S1", "S2",
                                       "EM", "N", "ORF3a6", "ORF7a7b8"))

# plot boxplot as inflammatory VS non-inf group
gCols <- c("Inflammatory Group" = "red", "Non-inflammatory Group" = "lightblue")

# CD4 Tcell stim
bPlot <- ggplot(data = ics,
                mapping = aes(x = stimulation,
                              y = pctpos_adj_CD8,
                              fill = iGroup)) +
                geom_boxplot(outlier.shape = NA, size = 0.3) +
            geom_point(position = position_jitterdodge(),
                        aes(fill = factor(iGroup)),
                        shape = 21,
                            width = 0.1,
                        size = 1) +
          scale_fill_manual(values = gCols) +
          labs(c = NULL, y = "%CD8+ T-cells") +
          theme_bw() +
          theme(strip.text = element_text(size = 12),
                axis.text.x = element_text(color = "black", size = 8),
                axis.ticks.x = element_blank(),
                axis.text.y = element_text(size = 7, color = "black"),
                axis.title.y = element_text(size = 15, color = "black"),
                axis.line = element_line(color = "grey"),
                panel.background = element_blank(),
                panel.grid = element_blank())
pdf(file = "Figure S3D.pdf", width = 10, height = 5, useDingbats = F)
print(bPlot)
dev.off()


# CD4 Tcell stim
bPlot <- ggplot(data = ics,
                mapping = aes(x = stimulation,
                              y = pctpos_adj_CD4,
                              fill = iGroup)) +
                geom_boxplot(outlier.shape = NA, size = 0.3) +
            geom_point(position = position_jitterdodge(),
                        aes(fill = factor(iGroup)),
                        shape = 21,
                            width = 0.1,
                        size = 1) +
          scale_fill_manual(values = gCols) +
          labs(c = NULL, y = "%CD4+ T-cells") +
          theme_bw() +
          theme(strip.text = element_text(size = 12),
                axis.text.x = element_text(color = "black", size = 8),
                axis.ticks.x = element_blank(),
                axis.text.y = element_text(size = 7, color = "black"),
                axis.title.y = element_text(size = 15, color = "black"),
                axis.line = element_line(color = "grey"),
                panel.background = element_blank(),
                panel.grid = element_blank())
pdf(file = "Figure S3C.pdf", width = 10, height = 5, useDingbats = F)
print(bPlot)
dev.off()

write.xlsx(ics, "FigurePDFs_Tables_Manuscript/Source data/Source Figure S3C S3D.xlsx", col.names = TRUE, row.names = FALSE, append = FALSE)



### RBD IgG estimates
abestDF <- read.csv("RBD IgG estimates.csv") %>%
           rename(PTID = PUBID) %>%
         filter(PTID %in% pivotDF$PTID) %>%
         mutate(clusterID = pivotDF$Cluster[match(PTID,
                           table = pivotDF$PTID)],
                Group = pivotDF$Group[match(PTID,
                               table = pivotDF$PTID)],
                iGroup = ifelse(clusterID %in% c("C4", "C5"),
                                    "Inflammatory Group",
                                    "Non-inflammatory Group")) %>%
        filter(!Group == "Uninfected") %>%
        select(-`COVID.PTID`)
# wilcox tests
wilcox.test(log10(IgGRBDd60_origscale) ~ iGroup, data = abestDF)

# plot boxplot
gCols <- c("Infected Recovered" = "black", "Infected PASC" = "magenta")
abestDF$Group <- factor(abestDF$Group, levels = names(gCols))

bPlot <- ggplot(data = abestDF,
                mapping = aes(x = clusterID,
                              y = log10(IgGRBDd60_origscale))) +
                geom_boxplot(outlier.shape = NA, size = 0.3) +
            geom_point(position = position_jitterdodge(),
                        aes(color = factor(Group)),
                            width = 0.1,
                        size = 1) +
          scale_color_manual(values = gCols) +
          labs(c = NULL, y = "log10 (RBD IgG titer estimates)",
               color = "Group") +
          theme_classic() +
          theme(strip.text = element_text(size = 10),
                axis.text.x = element_text(color = "black", size = 12),
                axis.ticks.x = element_blank(),
                axis.text.y = element_text(size = 7, color = "black"),
                axis.title.y = element_text(size = 15, color = "black"),
                axis.line = element_line(color = "black"))
pdf(file = "Fig 1B.pdf", width = 10, height = 5, useDingbats = F)
print(bPlot)
dev.off()



### Severity score difference between groups (ONLY AMONG PASC)
severityDF <- read_excel("Clinical activity score.xlsx", sheet = 1) %>%
              as.data.frame() %>%
              mutate(PUBID = pubids$PUBID[match(PTID, table = pubids$PTID)]) %>%
              filter(PUBID %in% pivotDF$PTID) %>%
              select(-Sample, -PTID) %>%
              mutate(Group = pivotDF$Group[match(PUBID,
                                                 table = pivotDF$PTID)],
                    clusterID = pivotDF$Cluster[match(PUBID,
                                                    table = pivotDF$PTID)],
                    iGroup = ifelse(clusterID %in% c("C4", "C5"),
                                                     "Inflammatory Group",
                                                     "Non-inflammatory Group")) %>%
              filter(Group == "Infected PASC")


# plot boxplot
severityDF$iGroup <- factor(severityDF$iGroup,
                            levels = c("Non-inflammatory Group",
                                       "Inflammatory Group"))

# box and jitter plot
gCols <- c("Uninfected" = "darkgrey",
          "Infected Recovered" = "black",
          "Infected PASC" = "magenta")

bPlot <- ggplot(data = severityDF,
                mapping = aes(x = clusterID,
                              y = `Clinical activity score`)) +
                geom_boxplot(outlier.shape = NA, size = 0.3) +
            geom_point(position = position_jitterdodge(),
                        aes(color = factor(Group)),
                        #shape = 21,
                            width = 0.1,
                        size = 1) +
          scale_color_manual(values = gCols) +
          theme_bw() +
          theme(strip.text = element_text(size = 10),
                axis.text.x = element_text(color = "black", size = 6),
                axis.ticks.x = element_blank(),
                axis.text.y = element_text(size = 7, color = "black"),
                axis.title.y = element_text(size = 15, color = "black"),
                axis.line = element_line(color = "grey"),
                panel.background = element_blank(),
                panel.grid = element_blank())
pdf(file = "Figure 1C.pdf", width = 6, height = 6, useDingbats = F)
print(bPlot)
dev.off()





### assessing BMI across clusters
plotDF <- tableS1 %>%
          select(PTID, Cluster, `BMI at enrollment`, Age, Group) %>%
          unique() %>%
          mutate(Group2 = ifelse(Cluster %in% c(4,5) &
                                 Group == "Infected PASC",
                                 "Inflammatory PASC",
                          ifelse(Cluster %in% c(2,3,4) &
                                 Group == "Infected PASC",
                                 "Non-inflammatory PASC",
                          Group)),
                iGroup = ifelse(Cluster %in% c(4,5),
                         "inflammatory", "non-inflammatory"),
                `BMI at enrollment` = as.numeric(`BMI at enrollment`),
                Age = as.numeric(Age))
gcols1 <- c("Uninfected" = "grey",
                "Infected Recovered" = "black",
                "Non-inflammatory PASC" = "blue",
                "Inflammatory PASC" = "red")

gcols2 <- c("Uninfected" = "grey",
           "Infected Recovered" = "black",
           "Infected PASC" = "magenta")

plotDF$Group <- factor(plotDF$Group, level = names(gcols2))
plotDF$Group2 <- factor(plotDF$Group2, level = names(gcols1))


bPlot <- ggplot(data = plotDF,
                mapping = aes(x = as.factor(Group),
                              y = `BMI at enrollment`)) +
          geom_boxplot(outlier.shape = NA,
                       size = 0.3) +
          geom_point(position = position_jitterdodge(jitter.width = 0.4),
                     aes(color = factor(Group2)),
                         size = 2) +
          scale_color_manual(values = gcols1) +
          labs(x = NULL, y = "BMI at enrollment",
               color = "Group") +
          theme_classic() +
          theme(axis.text.x = element_text(color = "black", size = 9),
                axis.ticks.x = element_blank(),
                axis.text.y = element_text(size = 15, color = "black"),
                panel.background = element_blank(),
                panel.grid = element_blank())
pdf(file = "Figure S4B.pdf", width = 6, height = 4, useDingbats = F)
print(bPlot)
dev.off()


bPlot <- ggplot(data = plotDF,
                mapping = aes(x = as.factor(Cluster),
                              y = `BMI at enrollment`)) +
          geom_boxplot(outlier.shape = NA,
                       size = 0.3) +
          geom_point(position = position_jitterdodge(jitter.width = 0.4),
                     aes(color = factor(Group)),
                         size = 2) +
          scale_color_manual(values = gcols2) +
          labs(x = "ClusterID",
               y = "BMI at enrollment",
               color = "Group") +
          theme_classic() +
          theme(axis.text.x = element_text(color = "black", size = 15),
                axis.ticks.x = element_blank(),
                axis.text.y = element_text(size = 15, color = "black"),
                panel.background = element_blank(),
                panel.grid = element_blank())
pdf(file = "Figure 1D.pdf", width = 6, height = 4, useDingbats = F)
print(bPlot)
dev.off()


### assessing Age across clusters
bPlot <- ggplot(data = plotDF,
                mapping = aes(x = as.factor(Cluster), y = Age)) +
          geom_boxplot(outlier.shape = NA,
                       size = 0.3) +
          geom_point(position = position_jitterdodge(jitter.width = 0.4),
                     aes(color = factor(Group)),
                     size = 2) +
          scale_color_manual(values = gcols2) +
          labs(x = "ClusterID",
               y = "Age",
               color = "Group") +
          theme_classic() +
          theme(strip.text = element_text(size = 8),
                axis.text.x = element_text(color = "black", size = 9),
                axis.ticks.x = element_blank(),
                axis.text.y = element_text(size = 12, color = "black"),
                axis.line = element_line(color = "grey"),
                panel.background = element_blank(),
                panel.grid = element_blank())
pdf(file = "Figure 1G.pdf", width = 6, height = 4, useDingbats = F)
print(bPlot)
dev.off()


bPlot <- ggplot(data = plotDF,
                mapping = aes(x = as.factor(Group), y = Age)) +
          geom_boxplot(outlier.shape = NA,
                       size = 0.3) +
          geom_point(position = position_jitterdodge(jitter.width = 0.4),
                     aes(color = factor(Group2)),
                         size = 2) +
          scale_color_manual(values = gcols1) +
          labs(x = NULL, y = "Age", color = "Group") +
          theme_classic() +
          theme(axis.text.x = element_text(color = "black", size = 9),
                axis.ticks.x = element_blank(),
                axis.text.y = element_text(size = 15, color = "black"),
                panel.background = element_blank(),
                panel.grid = element_blank())
pdf(file = "Figure S4A.pdf", width = 6, height = 4, useDingbats = F)
print(bPlot)
dev.off()



### check overlap with BMI and age related proteins in exisiting literature
## plot Venn
obesitymarkers <- read_excel("Obesity and age markers.xlsx", sheet = 1) %>%
                  as.data.frame() %>%
                 .$Protein %>%
                  unique()
agemarkers <- read_excel("Obesity and age markers.xlsx", sheet = 2) %>%
              as.data.frame() %>%
             .$Protein %>%
              unique()

ovLS <- list(obesitymarkers, agemarkers, clusterproteins)
names(ovLS) <- c("Known high BMI markers",
                 "Known high age markers",
                 "Inflammatory PASC markers")
                
# plot Venn diagram
library("VennDiagram")
library("SuperExactTest")

fileName = "BMI_Age_infPASC_venn.png"
venn <- venn.diagram(ovLS,
                     filename = fileName,
                     main = gsub(".png", "", fileName))


SeT <- supertest(ovLS, n = length(unique(unlist(ovLS))))
summSeT <- summary(SeT)
intersectionDF <- summSeT$Table %>% rownames_to_column(var = "overlapset")
write.table(intersectionDF,
            file = "Table S5.txt",
            sep = "\t",
            quote = F,
            row.names = F)



### assessing comorbidities
# stacked bar plot of each comorbiity proportions across clusters
tableS1 <- tableS1 %>%
           mutate(igroup = ifelse(Cluster %in% c(4,5),
                                  "Inflammatory group",
                                  "Non-inflammatory group"))
plotDF1 <- tableS1[, c(1, 8, 21:35)] %>% unique()
np_percluster <- as.data.frame(table(plotDF1$Cluster))

pdf(file = "Figure S4D.pdf", width = 4, height = 4)
for(i in c(3:(ncol(plotDF1)-1))) {

        plotDF2 <- plotDF1[, c(1,2,i)]
        
        # fisher test between inflammatory and non-inf groups
        print(colnames(plotDF1)[i])
        if(length(unique(plotDF2[, 3])) > 1) {
            print(fisher.test(table(plotDF2$igroup, plotDF2[, 3])))
                }
        
        
        tabDF <- as.data.frame(table(plotDF2$Cluster, plotDF2[, 3])) %>%
                 mutate(np = np_percluster$Freq[match(Var1,
                                          table = np_percluster$Var1)],
                        prop = (Freq/np)*100)

         # stacked bar of proportions
         fCols <- c("Yes" = "red", "No" = "black")
         plotbar <- ggplot(data = tabDF,
                            mapping =  aes(x = as.factor(Var1),
                                           y = prop,
                                           fill = Var2,
                                           label = round(prop, 1))) +
                     geom_bar(position = "stack", stat = "identity") +
                     geom_text(size = 3.5,
                               position = position_stack(vjust = 0.5),
                               color = "white") +
                     scale_fill_manual(values = fCols) +
                     labs(x = "Cluster number",
                          y = "%Participants",
                          fill = colnames(plotDF)[i]) +
                     theme_classic() +
                     theme(axis.text.x = element_text(size = 12, color = "black"),
                           axis.text.y = element_text(size = 12, color = "black"),
                           axis.line = element_line(color = "grey"),
                           panel.background = element_blank(),
                           panel.grid = element_blank())
         print(plotbar)
       }
dev.off()





#### read olink data
fullOlink <- read_excel("Table S2.xlsx", sheet = 1) %>% as.data.frame() %>%
             column_to_rownames(var = "Protein")

# subset full olink data on rule-in timepoints
Olink_rulein <- fullOlink[, pivotDF$PTID_visit, drop = F]



######## Identify markers of each Group: Uninfected, Recovered, PASC
groupIDs <- unique(pivotDF$Group)

clLS <- lapply(groupIDs, function(g) {
           print(g)
                                        
           # extract samples that belong to the cluster
            p1samples <- pivotDF %>% filter(Group == g) %>% .$PTID_visit
                                              
            # extract samples that belong all other clusters
            p2samples <- setdiff(pivotDF$PTID_visit, p1samples)
            
          # calculate delta of expression
          # difference in each protein expression between the selected cluster samples and all other samples
          olinkDF <- t(Olink_rulein)
                
          data.1 <- apply(X = olinkDF[p1samples, , drop = F],
                          MARGIN = 2,
                          FUN = function(x) mean(x))
          data.2 <- apply(X = olinkDF[p2samples, , drop = F],
                          MARGIN = 2,
                          FUN = function(x) mean(x))
          delta <- data.1-data.2
             
         # perform wilcox test: comparing protein expression in one cluster compared to all the other clusters
          group <- ifelse(rownames(olinkDF) %in% p1samples, "G1", "G2")
          group <- factor(group)
          p_val <- sapply(X = 1:ncol(x = olinkDF), FUN = function(x) {
                      return(wilcox.test(as.numeric(olinkDF[, x]) ~
                                         group)$p.value)
                                      })
          to.return <- data.frame(p_val = p_val,
                                  protein = colnames(olinkDF)) %>%
                       mutate(mCluster = as.numeric(data.1)
                                         [match(protein, names(data.1))],
                              mOtherClusters = as.numeric(data.2)
                                             [match(protein,names(data.2))],
                              delta = as.numeric(delta)
                                      [match(protein, names(delta))],
                              Group = g,
                              adjp_val = p.adjust(p_val, method = "BH")) %>%
                       filter(delta > 0 & p_val < 0.05)
          return(value = to.return)
        })
clDF_groups <- do.call(rbind, clLS) %>% arrange(Group, p_val)
        
# write table
clDF_groups <- clDF_groups %>% arrange(Group, p_val)
fileName <- "Table S3.txt"
write.table(clDF_groups,
            file = fileName,
            row.names = FALSE,
            sep = "\t",
            quote = FALSE)
    


## heatmap of significant DEPs per group
# take top 50 deps per cluster
topdeps <- clDF_groups %>%
           group_by(Group) %>%
           top_n(50, desc(p_val)) %>%
           as.data.frame() %>%
           .$protein
clDF_top <- clDF_groups %>%
            filter(protein %in% topdeps) %>%
            arrange(Group, p_val)

# subset on signifincat proteins slected above
subOlink <- Olink_rulein[topdeps, , drop = F]
subOlink <- t(scale(t(subOlink))) %>% as.data.frame()

# heatmap row annotation
matannot_row <- clDF_top %>%
                select(protein, Group) %>%
                column_to_rownames(var = "protein")

# heatmap column annotation
matannot_col <- pivotDF[, c(1, 7)] %>%
                column_to_rownames(var = "PTID_visit")
matannot_col <- matannot_col[colnames(subOlink), , drop = F]
table(colnames(subOlink) == rownames(matannot_col))

# column order
corder <- matannot_col %>%
          rownames_to_column() %>%
          arrange(Group) %>%
          .$rowname

# set limits
limit = 3

# annotation colors
ann_colors <- list(Group = c("Uninfected" = "grey",
                             "Infected Recovered" = "black",
                              "Infected PASC" = "magenta"))

# define custom color palatte
colorPalette <- c("purple", "black", "yellow")
colorPalette <- colorRampPalette(colors = colorPalette)(100)

# plot
outFile <- "Figure S2.pdf"
pdf(file = outFile, width = 10, height = 10)
        
pheatmap(mat = subOlink[rownames(matannot_row), corder],
        color = colorPalette,
        breaks = c(min(subOlink),
                   seq(from = -1*limit,
                       to = limit,
                       length.out = 99),
                   max(subOlink)),
        cellwidth = 2,
        cellheight = 5,
        cluster_cols = FALSE,
        cluster_rows = FALSE,
        annotation_col = matannot_col,
        annotation_row = matannot_row,
        show_rownames = TRUE,
        show_colnames = FALSE,
        annotation_colors = ann_colors,
        border_color = NA,
        fontsize_row = 5,
        fontsize = 7)
dev.off()




### Perform ssGSEA on modules and rank modules by top expressed in each cluster
# read modules that define clustering
sigmods <- read_excel("c2-Clustering-FTP-Longitudinal.xlsx", sheet = 3) %>%
           as.data.frame() %>%
           mutate(name = tolower(name), name = gsub("_" ," ", name))
gsLS <- sigmods %>% .$genes %>% strsplit(",")
names(gsLS) <- sigmods$module


# run ssGSEA: module expression across full data set
set.seed(20220122)
ssgsea_res <- gsva(data.matrix(fullOlink),
                   gsLS,
                   method = "ssgsea",
                   abs.ranking = FALSE,
                   min.sz = 2,
                   max.sz = 2000,
                   ssgsea.norm = T)
ssgsea_res <- t(scale(t(ssgsea_res)))


## Identify markers (modules) of each sub-cluster
clusterIDs <- unique(pivotDF$Cluster)

clLS <- lapply(clusterIDs, function(n) {
           print(n)
                                        
           # extract samples that belong to the cluster
            p1samples <- pivotDF %>% filter(Cluster == n) %>% .$PTID_visit
                                              
            # extract samples that belong all other clusters
            p2samples <- setdiff(pivotDF$PTID_visit, p1samples)
            
          # calculate delta of expression
          # difference in each protein expression between the selected cluster samples and all other samples
          ssgseaDF <- t(ssgsea_res)
                
          data.1 <- apply(X = ssgseaDF[p1samples, , drop = F],
                          MARGIN = 2,
                          FUN = function(x) mean(x))
          data.2 <- apply(X = ssgseaDF[p2samples, , drop = F],
                          MARGIN = 2,
                          FUN = function(x) mean(x))
          delta <- data.1-data.2
             
         # perform wilcox test: comparing protein expression in one cluster compared to all the other clusters
          group <- ifelse(rownames(ssgseaDF) %in% p1samples, "G1", "G2")
          group <- factor(group)
          p_val <- sapply(X = 1:ncol(x = ssgseaDF), FUN = function(x) {
                      return(wilcox.test(as.numeric(ssgseaDF[, x]) ~
                                         group)$p.value)
                                      })
          to.return <- data.frame(p_val = p_val,
                                  module = colnames(ssgseaDF)) %>%
                       mutate(mCluster = as.numeric(data.1)
                                         [match(module, names(data.1))],
                              mOtherClusters = as.numeric(data.2)
                                             [match(module,names(data.2))],
                              delta = as.numeric(delta)
                                      [match(module, names(delta))],
                              clusterID = n,
                              adjp_val = p.adjust(p_val, method = "BH")) %>%
                       filter(delta > 0 & p_val < 0.05)
          return(value = to.return)
        })
clDF1 <- do.call(rbind, clLS) %>%
         arrange(clusterID, p_val) %>%
         mutate(gsName = sigmods$name[match(module, table = sigmods$module)],
                gsName = tolower(gsName),
                gsName = gsub("_" ," ", gsName))
                 


######## Identify individual DEPs of each sub-cluster
clusterIDs <- unique(pivotDF$Cluster)

clLS <- lapply(clusterIDs, function(n) {
           print(n)
                                        
           # extract samples that belong to the cluster
            p1samples <- pivotDF %>% filter(Cluster == n) %>% .$PTID_visit
                                              
            # extract samples that belong all other clusters
            p2samples <- setdiff(pivotDF$PTID_visit, p1samples)
            
          # calculate delta of expression
          # difference in each protein expression between the selected cluster samples and all other samples
          olinkDF <- t(Olink_rulein)
                
          data.1 <- apply(X = olinkDF[p1samples, , drop = F],
                          MARGIN = 2,
                          FUN = function(x) mean(x))
          data.2 <- apply(X = olinkDF[p2samples, , drop = F],
                          MARGIN = 2,
                          FUN = function(x) mean(x))
          delta <- data.1-data.2
             
         # perform wilcox test: comparing protein expression in one cluster compared to all the other clusters
          group <- ifelse(rownames(olinkDF) %in% p1samples, "G1", "G2")
          group <- factor(group)
          p_val <- sapply(X = 1:ncol(x = olinkDF), FUN = function(x) {
                      return(wilcox.test(as.numeric(olinkDF[, x]) ~
                                         group)$p.value)
                                      })
          to.return <- data.frame(p_val = p_val,
                                  protein = colnames(olinkDF)) %>%
                       mutate(mCluster = as.numeric(data.1)
                                         [match(protein, names(data.1))],
                              mOtherClusters = as.numeric(data.2)
                                             [match(protein,names(data.2))],
                              delta = as.numeric(delta)
                                      [match(protein, names(delta))],
                              clusterID = n,
                              adjp_val = p.adjust(p_val, method = "BH")) %>%
                       filter(delta > 0 & p_val < 0.05)
          return(value = to.return)
        })
clDF2 <- do.call(rbind, clLS) %>% arrange(clusterID, p_val)


### correlation between individual olink markers and BMI
clusterproteins <- clDF2 %>%
                   filter(adjp_val < 0.05) %>%
                   filter(clusterID %in% c("H4", "H5")) %>%
                  .$protein %>%
                   unique()

bmiDF <- tableS1 %>%
         select(PTID_visit, Group, Cluster,`BMI at enrollment`) %>%
         filter(PTID_visit %in% colnames(Olink_rulein)) %>%
         filter(!`BMI at enrollment` %in% NA) %>%
         filter(!Group == "Uninfected") %>%
         mutate(`BMI at enrollment` = as.numeric(`BMI at enrollment`)) %>%
         column_to_rownames(var = "PTID_visit")
subDF <- Olink_rulein[clusterproteins, rownames(bmiDF), drop = F]
table(rownames(bmiDF) == colnames(subDF))


subDF2 <- subDF %>%
         rownames_to_column(var = "protein") %>%
         gather(PTID_visit, NPX, -protein) %>%
         mutate(BMI = bmiDF$"BMI at enrollment"[match(PTID_visit,
                            table = rownames(bmiDF))],
                Group = bmiDF$Group[match(PTID_visit,
                            table = rownames(bmiDF))])
corDF <- subDF2 %>%
         group_by(protein) %>%
         do(p = cor.test(.$NPX, .$BMI, method = "spearman")$p.value,
            rho = cor.test(.$NPX, .$BMI, method = "spearman")$estimate) %>%
         mutate(p = unlist(p), rho = unlist(rho)) %>%
         mutate(padj = p.adjust(p, method = "BH")) %>%
         as.data.frame() %>%
         filter(padj < 0.05 & rho > 0.49) %>%
         arrange(padj)

# heatmap of significant correlates
library("ComplexHeatmap")
library("circlize")
subDF2_sig <- subDF[corDF$protein, , drop = F]
subDF2_sig <- t(scale(t(subDF2_sig)))
bmiDF2 <- bmiDF %>%
          rownames_to_column() %>%
          arrange(`BMI at enrollment`) %>%
          column_to_rownames(var = "rowname")
subDF2_sig <- subDF2_sig[, rownames(bmiDF2), drop = F]
table(rownames(bmiDF2) == colnames(subDF2_sig))
colannot <- HeatmapAnnotation(df = bmiDF2 %>%
                                   dplyr::select(Cluster),
                              BMI = anno_barplot(bmiDF2$"BMI at enrollment"),
                              col = list(Cluster = c(`1` = "grey",
                                                    `2` = "black",
                                                    `3` = "blue",
                                                    `4` = "red",
                                                    `5` = "orange")),
                              gp = gpar(col = "black"))

htmap <- Heatmap(subDF2_sig,
                 cluster_rows = F,
                 cluster_columns = F,
                 na_col = "grey",
                 col = colorRamp2(c(-3,0,3),
                                  c("purple",
                                    "black",
                                    "yellow")),
                 row_names_max_width = unit(10, "cm"),
                 show_column_names = F,
                 column_names_gp = gpar(fontsize = 5),
                 row_names_gp = gpar(fontsize = 10),
                 top_annotation = colannot,
                 heatmap_legend_param = list(title = "NPX(z-score)",
                                             heatmap_legend_side = "right"))
outFile <- "Figure 1E.pdf"
pdf(file = outFile, width = 5, height = 4)
htmap
dev.off()


write.xlsx(subDF2_sig, "FigurePDFs_Tables_Manuscript/Source data/Source Figure 1E.xlsx",
  col.names = TRUE, row.names = TRUE, append = FALSE)



### correlation between individual modules and BMI
modulesDF1 <- ssgsea_res
modulesDF1 <- modulesDF1[,  rownames(bmiDF2), drop = F]
modulesDF2 <- modulesDF1 %>%
              as.data.frame() %>%
              rownames_to_column(var = "module") %>%
              gather(PTID_visit, ssGSEA, -module) %>%
              mutate(BMI = bmiDF2$"BMI at enrollment"[match(PTID_visit,
                                table = rownames(bmiDF2))],
                clusterID = bmiDF2$Cluster[match(PTID_visit,
                            table = rownames(bmiDF2))])
corDF <- modulesDF2 %>%
         group_by(module) %>%
         do(p = cor.test(.$ssGSEA, .$BMI, method = "spearman")$p.value,
            rho = cor.test(.$ssGSEA, .$BMI, method = "spearman")$estimate) %>%
         mutate(p = unlist(p), rho = unlist(rho)) %>%
         mutate(padj = p.adjust(p, method = "BH")) %>%
         as.data.frame() %>%
         filter(padj < 0.05 & rho > 0) %>%
         arrange(padj)

# correlation scatter plots of top significant correlation
plotDF <- modulesDF2 %>% filter(module == "module6")
ccols <- c(`1` = "grey",
            `2` = "black",
            `3` = "blue",
            `4` = "red",
            `5` = "orange")
cplot <- ggplot(data = plotDF,
                 mapping = aes(x = ssGSEA, y = BMI,
                               color = as.factor(clusterID))) +
         geom_point(size = 2) +
         scale_color_manual(values = ccols) +
         labs(x = "ssGSEA score: Leptin signaling pathway (M6)",
              y = "BMI at enrollment") +
         geom_smooth(method = "lm", color = "black", se = TRUE) +
         theme_classic() +
         theme(axis.text.x = element_text(size = 10, color = "black"),
               axis.text.y = element_text(size = 10, color = "black"))
pdf(file = "Figure 1F.pdf", width = 4, height = 4)
print(cplot)
dev.off()


##### correlation between individual proteins and age
ageDF <- tableS1 %>%
         select(PTID_visit, Group, Age, Cluster) %>%
         filter(PTID_visit %in% colnames(Olink_rulein)) %>%
         filter(!Group == "Uninfected") %>%
         arrange(Age) %>%
         column_to_rownames(var = "PTID_visit")
subDF <- Olink_rulein[clusterproteins, rownames(ageDF), drop = F]
table(rownames(ageDF) == colnames(subDF))

subDF2 <- subDF %>%
         rownames_to_column(var = "protein") %>%
         gather(PTID_visit, NPX, -protein) %>%
         mutate(Age = ageDF$Age[match(PTID_visit,table = rownames(ageDF))],
                Group = ageDF$Group[match(PTID_visit,
                            table = rownames(ageDF))])
corDF <- subDF2 %>%
         group_by(protein) %>%
         do(p = cor.test(.$NPX, .$Age, method = "spearman")$p.value,
            rho = cor.test(.$NPX, .$Age, method = "spearman")$estimate) %>%
         mutate(p = unlist(p), rho = unlist(rho)) %>%
         mutate(padj = p.adjust(p, method = "BH")) %>%
         as.data.frame() %>%
         filter(padj < 0.05 & rho >= 0.45) %>%
         arrange(padj)

# heatmap of significant correlates
subDF2_sig <- subDF[corDF$protein, , drop = F]
subDF2_sig <- t(scale(t(subDF2_sig)))
ageDF2 <- ageDF %>%
          rownames_to_column() %>%
          arrange(Age) %>%
          column_to_rownames(var = "rowname")
subDF2_sig <- subDF2_sig[, rownames(ageDF2), drop = F]
table(rownames(ageDF2) == colnames(subDF2_sig))
colannot <- HeatmapAnnotation(df = ageDF2 %>%
                                    dplyr::select(Cluster),
                              Age = anno_barplot(ageDF2$Age),
                              col = list(`Rule.in.cluster` = c(`1` = "grey",
                                                            `2` = "black",
                                                            `3` = "blue",
                                                            `4` = "red",
                                                            `5` = "orange")),
                              gp = gpar(col = "black"))

htmap <- Heatmap(subDF2_sig,
                 cluster_rows = F,
                 cluster_columns = F,
                 na_col = "grey",
                 col = colorRamp2(c(-3,0,3),
                                  c("purple",
                                    "black",
                                    "yellow")),
                 row_names_max_width = unit(10, "cm"),
                 show_column_names = F,
                 column_names_gp = gpar(fontsize = 5),
                 row_names_gp = gpar(fontsize = 5),
                 top_annotation = colannot,
                 heatmap_legend_param = list(title = "NPX(z-score)",
                                             heatmap_legend_side = "right"))
outFile <- "Figure 1H.pdf"
pdf(file = outFile, width = 5, height = 4)
htmap
dev.off()



### correlation between individual modules and age
modulesDF1 <- ssgsea_res
modulesDF1 <- modulesDF1[,  rownames(ageDF), drop = F]
modulesDF2 <- modulesDF1 %>%
              as.data.frame() %>%
              rownames_to_column(var = "module") %>%
              gather(PTID_visit, ssGSEA, -module) %>%
              mutate(Age = ageDF$Age[match(PTID_visit,
                                table = rownames(ageDF))],
                     clusterID = ageDF$Cluster[match(PTID_visit,
                                    table = rownames(ageDF))])
corDF <- modulesDF2 %>%
         group_by(module) %>%
         do(p = cor.test(.$ssGSEA, .$Age, method = "spearman")$p.value,
            rho = cor.test(.$ssGSEA, .$Age, method = "spearman")$estimate) %>%
         mutate(p = unlist(p), rho = unlist(rho)) %>%
         mutate(padj = p.adjust(p, method = "BH")) %>%
         as.data.frame() %>%
         filter(padj < 0.05 & rho > 0) %>%
         arrange(padj)

# correlation scatter plots of top significant correlation
plotDF <- modulesDF2 %>% filter(module == "module32")
cplot <- ggplot(data = plotDF,
                 mapping = aes(x = ssGSEA, y = Age,
                                color = as.factor(clusterID))) +
         geom_point(size = 2) +
         scale_color_manual(values = ccols) +
         labs(x = "ssGSEA score: Type II interferon signaling IFNG (M32)",
              y = "Age at enrollment") +
         geom_smooth(method = "lm", color = "black", se = TRUE) +
         theme_classic() +
         theme(axis.text.x = element_text(size = 10, color = "black"),
               axis.text.y = element_text(size = 10, color = "black"))
pdf(file = "Figure 1I.pdf", width = 4, height = 4)
print(cplot)
dev.off()



##### network of significant modules that sigificantly define groups
graphDF <- clDF1 %>%
           filter(adjp_val < 0.05) %>%
           mutate(lp = -log10(adjp_val),
                  clusterID = gsub("C", "Cluster", clusterID)) %>%
           select(clusterID, gsName, lp) %>%
           rename(from = clusterID, to = gsName) %>%
           transform(count = table(to)[to]) %>%
           mutate(`count.Freq` = (100/`count.Freq`)/100)

# create igraph object
g <- graph.data.frame(graphDF, directed = FALSE)
g <- simplify(g, remove.loops = TRUE)

# create vertices as pies
vpieDF <- rbind(data.frame(`count.Freq` = rep(1, 4),
                           `count.to` = V(g)$name[1:4]),
                graphDF[, c("count.to", "count.Freq")])
values <- unstack(vpieDF)
values <- values[V(g)$name]

# generate vertex colors based on p-values
minVal <- min(graphDF$lp)
maxVal <- max(graphDF$lp)

colorPalette <- breaks <- NULL
if (0 > minVal & 0 < maxVal) {
    breakVal <- min(abs(range(graphDF$lp)))
    breaks <- c(floor(minVal),
                seq(from = -breakVal,
                    to = breakVal,
                    length.out = 99),
                ceiling(maxVal))
    colorPalette <- c("#D73027", "#FC8D59", "#FEE090", "#FFFFFF",
                      "#E0F3F8", "#91BFDB", "#4575B4")
    colorPalette <- colorRampPalette(colors = rev(colorPalette))(100)
    } else {
    if (0 <= minVal) {
        breaks <- seq(from = 0, to = maxVal, length.out = 101)
        colorPalette <- c("#D73027", "#FC8D59", "#FEE090","#FFFFFF")
        colorPalette <- colorRampPalette(colors = rev(colorPalette))(100)
    } else {
        breaks <- seq(from = minVal, to = 0, length.out = 101)
        colorPalette <- c("#FFFFFF", "#E0F3F8", "#91BFDB","#4575B4")
        colorPalette <- colorRampPalette(colors = rev(colorPalette))(100)
    }
 }
keyColor <- findInterval(graphDF$lp, vec = breaks, all.inside = TRUE)
keyColor <- colorPalette[keyColor]

# create vertex pie colors object
vpieDF$vcol <- c(rep("grey", 4), keyColor)
values_col <- unstack(vpieDF[, c("vcol", "count.to")])
values_col <- values_col[V(g)$name]

# tkplot graph to get custom coordinates
# for arranging nodes
graphDF %>%
filter(`count.Freq` == 0.5) %>%
group_by(to) %>%
summarize(m = median(lp)) %>%
as.data.frame() %>%
arrange(desc(m))

graphDF %>%
filter(`count.Freq` == 1) %>%
group_by(from, to) %>%
summarize(m = median(lp)) %>%
as.data.frame() %>%
arrange(from, desc(m))

tkplot(g, vertex.color = "pink", vertex.size = 8)

# get node coordinates generated by tkplot
posMat <- tkplot.getcoords(tkp.id = 1)

# save graph
outFile <- "Figure 2A.pdf"
pdf(file = outFile, width = 20, height = 20)
plot(g,
     vertex.shape = "pie",
     vertex.pie = values,
     vertex.pie.color = values_col,
     vertex.size = 6,
     vertex.label.color = "black",
     vertex.frame.color = "black",
     vertex.label.cex = 0.8,
     layout = posMat)
dev.off()




#### Box and jitter plot of cluster biomarkers modules
selModules <- unique(clDF1$module)

plotDF <- ssgsea_res[selModules, pivotDF$PTID_visit, drop = F] %>%
          as.data.frame() %>%
          rownames_to_column(var = "module") %>%
          gather(PTID_visit, ssgsea, -module) %>%
          mutate(clusterID = pivotDF$Cluster[match(PTID_visit,
                            table = pivotDF$PTID_visit)],
                 Group = pivotDF$Group[match(PTID_visit,
                                        table = pivotDF$PTID_visit)]) %>%
            mutate(gsName = sigmods$name[match(module,
                                            table = sigmods$module)])

# box and jitter plot
gCols <- c("Uninfected" = "darkgrey",
          "Infected Recovered" = "black",
          "Infected PASC" = "magenta")
plotDF$Group <- factor(plotDF$Group, levels = names(gCols))

bPlot <- ggplot(data = plotDF,
                mapping = aes(x = clusterID, y = ssgsea)) +
          geom_boxplot(outlier.shape = NA,
                       size = 0.3) +
          geom_point(position = position_jitterdodge(),
                      aes(color = factor(Group)),
                      width = 0.1,
                      size = 1) +
          scale_color_manual(values = gCols) +
          facet_wrap(~gsName, scales = "free_y", ncol = 7) +
          labs(x = NULL,
               y = "ssGSEA score",
               color = "Group") +
          theme_bw() +
          theme(strip.text = element_text(size = 8),
                axis.text.x = element_text(color = "black", size = 9),
                axis.ticks.x = element_blank(),
                axis.text.y = element_text(size = 12, color = "black"),
                axis.line = element_line(color = "grey"),
                panel.background = element_blank(),
                panel.grid = element_blank())
pdf(file = "Figure S5.pdf", width = 35, height = 23, useDingbats = F)
print(bPlot)
dev.off()



##### Correlellogram on top modules defning C4 and C5
selmodules <- c("wp type ii interferon signaling ifng",
                "pid il27 pathway",
                "pid il1 pathway",
                "reactome regulation of ifna signaling",
                "reactome tnf signaling",
                "biocarta tid pathway",
                "wp il18 signaling pathway",
                "pid nfkappab canonical pathway")

subs_sgsea_res <- ssgsea_res
rownames(subs_sgsea_res) <- sigmods$name[match(rownames(subs_sgsea_res),
                                                table = sigmods$module)]
subs_sgsea_res <- subs_sgsea_res[selmodules, pivotDF$PTID_visit, drop = F]

# correlation tests
cor.test(as.numeric(subslea["pid nfkappab canonical pathway", ]),
         as.numeric(subslea["wp type ii interferon signaling ifng", ]),
         method = "spearman")

corDF <- cor(t(subs_sgsea_res), method = "spearman")
corDF <- round(corDF, 2)

# Get lower triangle of the correlation matrix
get_lower_tri <- function(cormat){
                    cormat[upper.tri(cormat)] <- NA
                    return(cormat)
                  }
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
                    cormat[lower.tri(cormat)]<- NA
                    return(cormat)
                  }

# reordering the correlation matrix
reorder_cormat <- function(cormat){
                    # Use correlation between variables as distance
                     dd <- as.dist((1-cormat)/2)
                     hc <- hclust(dd)
                     cormat <- cormat[hc$order, hc$order]
                }
corDF <- reorder_cormat(corDF)
upper_tri <- get_upper_tri(corDF)

# Melt the correlation matrix
library(reshape2)
melted_cormat <- melt(upper_tri, na.rm = TRUE)

# correlation Heatmap
library(viridis)
cplot <- ggplot(data = melted_cormat,
                aes(Var2, Var1, fill = value)) +
         geom_tile(color = "white")+
         scale_fill_gradient2(low = "blue",
                              high = "red",
                              mid = "white",
                              midpoint = 0,
                              limit = c(-1,1),
                              space = "Lab",
                              name = "Spearman\nCorrelation") +
         theme_classic() +
         theme(axis.text.x = element_text(size = 7,
                                          angle = 90,
                                          hjust = 1,
                                          color = "black"),
               axis.text.y = element_text(size = 7,
                                         color = "black")) +
         geom_text(aes(Var2, Var1, label = value),
                   color = "black",
                   size = 2)
pdf(file = "Figure S6.pdf", width = 5, height = 4)
print(cplot)
dev.off()



##### Network diagram of cluster specific proteins
# graph data frame to plot cytokines and cytokine receptors
# select cytokines and cytokine receptors
selMarkers <- rownames(fullOlink)[grep("^IL|EBI3|^CCL|^CXC|IFN|^TNF|^C2|^C4BPB|CSF1|CSF3|TGFB1$|TGFBR2|TGFBR3", rownames(fullOlink))]
selMarkers <- selMarkers[-grep("RN|ILKAP|C2CD2L", selMarkers)]

# make graph data frame
graphDF <- clDF2 %>%
           filter(adjp_val < 0.05) %>%
           filter(protein %in% selMarkers) %>%
           mutate(lp = -log10(p_val),
                  clusterID = gsub("C", "Cluster", clusterID)) %>%
           select(clusterID, protein, lp) %>%
           rename(from = clusterID, to = protein) %>%
           transform(count = table(to)[to]) %>%
           mutate(`count.Freq` = (100/`count.Freq`)/100)

# create igraph object
g <- graph.data.frame(graphDF, directed = FALSE)
g <- simplify(g, remove.loops = TRUE)

# create vertices as pies
vpieDF <- rbind(data.frame(`count.Freq` = rep(1, 2),
                           `count.to` = V(g)$name[1:2]),
                graphDF[, c("count.to", "count.Freq")])
values <- unstack(vpieDF)
values <- values[V(g)$name]

# generate vertex colors based on p-values
minVal <- min(graphDF$lp)
maxVal <- max(graphDF$lp)

colorPalette <- breaks <- NULL
if (0 > minVal & 0 < maxVal) {
    breakVal <- min(abs(range(graphDF$lp)))
    breaks <- c(floor(minVal),
                seq(from = -breakVal,
                    to = breakVal,
                    length.out = 99),
                ceiling(maxVal))
    colorPalette <- c("#D73027", "#FC8D59", "#FEE090", "#FFFFFF",
                      "#E0F3F8", "#91BFDB", "#4575B4")
    colorPalette <- colorRampPalette(colors = rev(colorPalette))(100)
    } else {
    if (0 <= minVal) {
        breaks <- seq(from = 0, to = maxVal, length.out = 101)
        colorPalette <- c("#D73027", "#FC8D59", "#FEE090","#FFFFFF")
        colorPalette <- colorRampPalette(colors = rev(colorPalette))(100)
    } else {
        breaks <- seq(from = minVal, to = 0, length.out = 101)
        colorPalette <- c("#FFFFFF", "#E0F3F8", "#91BFDB","#4575B4")
        colorPalette <- colorRampPalette(colors = rev(colorPalette))(100)
    }
 }
keyColor <- findInterval(graphDF$lp, vec = breaks, all.inside = TRUE)
keyColor <- colorPalette[keyColor]

# create vertex pie colors object
vpieDF$vcol <- c(rep("grey", 2), keyColor)
values_col <- unstack(vpieDF[, c("vcol", "count.to")])
values_col <- values_col[V(g)$name]

# tkplot graph to get custom coordinates
# for arranging nodes
graphDF %>%
filter(`count.Freq` == 0.5) %>%
group_by(to) %>%
summarize(m = median(lp)) %>%
as.data.frame() %>%
arrange(desc(m))

graphDF %>%
filter(`count.Freq` == 1) %>%
group_by(from, to) %>%
summarize(m = median(lp)) %>%
as.data.frame() %>%
arrange(from, desc(m))

tkplot(g, vertex.color = "pink", vertex.size = 9)

# get node coordinates generated by tkplot
posMat <- tkplot.getcoords(tkp.id = 12)

# save graph
outFile <- "Figure 3A.pdf"
pdf(file = outFile, width = 7, height = 7)
plot(g,
     vertex.shape = "pie",
     vertex.pie = values,
     vertex.pie.color = values_col,
     vertex.size = 10,
     vertex.label.color = "black",
     vertex.frame.color = "black",
     vertex.label.cex = 0.8,
     layout = posMat)
dev.off()


## Network of top 30 (full olink panel) proteins of each cluster
# graph data frame to plot top N deps (ranked by adjp_val) per cluster
graphDF <- clDF2 %>%
           filter(adjp_val < 0.05) %>%
           group_by(clusterID) %>%
           top_n(30, desc(adjp_val)) %>%
           as.data.frame() %>%
           mutate(lp = -log10(adjp_val),
                  clusterID = gsub("C", "Cluster", clusterID)) %>%
           select(clusterID, protein, lp) %>%
           rename(from = clusterID, to = protein) %>%
           transform(count = table(to)[to]) %>%
           mutate(`count.Freq` = (100/`count.Freq`)/100)

# create igraph object
g <- graph.data.frame(graphDF, directed = FALSE)
g <- simplify(g, remove.loops = TRUE)

# generate vertex colors based on p-values
minVal <- min(graphDF$lp)
maxVal <- max(graphDF$lp)

colorPalette <- breaks <- NULL
if (0 > minVal & 0 < maxVal) {
    breakVal <- min(abs(range(graphDF$lp)))
    breaks <- c(floor(minVal),
                seq(from = -breakVal,
                    to = breakVal,
                    length.out = 99),
                ceiling(maxVal))
    colorPalette <- c("#D73027", "#FC8D59", "#FEE090", "#FFFFFF",
                      "#E0F3F8", "#91BFDB", "#4575B4")
    colorPalette <- colorRampPalette(colors = rev(colorPalette))(100)
    } else {
    if (0 <= minVal) {
        breaks <- seq(from = 0, to = maxVal, length.out = 101)
        colorPalette <- c("#D73027", "#FC8D59", "#FEE090","#FFFFFF")
        colorPalette <- colorRampPalette(colors = rev(colorPalette))(100)
    } else {
        breaks <- seq(from = minVal, to = 0, length.out = 101)
        colorPalette <- c("#FFFFFF", "#E0F3F8", "#91BFDB","#4575B4")
        colorPalette <- colorRampPalette(colors = rev(colorPalette))(100)
    }
 }
keyColor <- findInterval(graphDF$lp, vec = breaks, all.inside = TRUE)
keyColor <- colorPalette[keyColor]
keyColor <- c(rep("grey", 4), keyColor)

# tkplot graph to get custom coordinates
# for arranging nodes
tkplot(g, vertex.color = "pink", vertex.size = 8)

# get node coordinates generated by tkplot
posMat <- tkplot.getcoords(tkp.id = 1)

# save graph
outFile <- "Figure S8.pdf"
pdf(file = outFile, width = 7, height = 7, useDingbats = F)
plot(g,
     vertex.color = keyColor,
     vertex.size = 7,
     vertex.label.color = "black",
     vertex.frame.color = "black",
     vertex.label.cex = 0.8,
     layout = posMat)
dev.off()




### Box and jitter plot of cluster biomarkers (top 30 per cluster ranked by adj-pval)
selproteins <- clDF2 %>%
               filter(adjp_val < 0.05) %>%
               group_by(clusterID) %>%
               top_n(30, desc(adjp_val)) %>%
              .$protein

plotDF <- Olink_rulein[selproteins, , drop = F] %>%
          rownames_to_column(var = "protein") %>%
          gather(PTID_visit, NPX, -protein) %>%
          mutate(clusterID = pivotDF$Cluster[match(PTID_visit,
                            table = pivotDF$PTID_visit)],
                 Group = pivotDF$Group[match(PTID_visit,
                                    table = pivotDF$PTID_visit)])

plotDF$protein <- factor(plotDF$protein, levels = selproteins)

# box and jitter plot
gCols <- c("Uninfected" = "darkgrey",
          "Infected Recovered" = "black",
          "Infected PASC" = "magenta")
plotDF$Group <- factor(plotDF$Group, levels = names(gCols))

bPlot <- ggplot(data = plotDF,
                mapping = aes(x = clusterID, y = NPX)) +
          geom_boxplot(outlier.shape = NA,
                       size = 0.3) +
          geom_point(position = position_jitterdodge(),
                     aes(color = factor(Group)),
                         #shape = 21,
                     width = 0.1,
                     size = 1.5) +
          #scale_x_discrete(limits = rev(corder)) +
          scale_color_manual(values = gCols) +
          facet_wrap(~protein, scales = "free", ncol = 6) +
          labs(x = NULL,
               y = "Normalized Protein Expression",
               color = "Group") +
          theme_bw() +
          theme(strip.text = element_text(size = 12),
                axis.text.x = element_text(color = "black", size = 9),
                axis.ticks.x = element_blank(),
                axis.text.y = element_text(size = 12, color = "black"),
                axis.line = element_line(color = "grey"),
                panel.background = element_blank(),
                panel.grid = element_blank())
pdf(file = "Figure S7.pdf", width = 35, height = 35, useDingbats = F)
print(bPlot)
dev.off()



## heatmap of top 50 significant DEPs per cluster ID (C4 and C5)
clDF_top <- clDF2 %>%
            filter(clusterID %in% c("H4", "H5")) %>%
            group_by(clusterID) %>%
            top_n(50, desc(adjp_val)) %>%
            as.data.frame()

# subset on signifincat proteins selected above
subOlink <- Olink_rulein[clDF_top$protein, , drop = F]
subOlink <- t(scale(t(subOlink))) %>% as.data.frame()


# arraging by proteins that are uniquely and commonly expressed across different sub-groups
dupdeps <- clDF_top$protein[duplicated(clDF_top$protein)]
commclusts <- clDF_top %>%
              group_by(protein) %>%
              summarize(cl = paste(clusterID, collapse = "_")) %>%
              as.data.frame()

# heatmap row annotation
matannot_row <- clDF_top %>%
                select(protein, clusterID) %>%
                mutate(clusterID2 = ifelse(protein %in% dupdeps,
                                            commclusts$cl[match(protein,
                                                table = commclusts$protein)],
                                           clusterID)) %>%
                select(protein, clusterID2) %>%
                unique() %>%
                arrange(clusterID2) %>%
                mutate(clusterID2 = as.character(clusterID2)) %>%
                column_to_rownames(var = "protein")

# heatmap column annotation
matannot_col <- pivotDF[, c(1, 4, 7, 10:17)] %>%
                column_to_rownames(var = "PTID_visit")
matannot_col <- matannot_col[colnames(subOlink), , drop = F]
table(colnames(subOlink) == rownames(matannot_col))

# column order
corder <- matannot_col %>%
          rownames_to_column() %>%
          arrange(Cluster, rowname) %>%
          .$rowname

# set limits
limit = 3

# annotation colors
ann_colors <- list(Cluster = c("C1" = "black",
                               "C2" = "blue",
                               "C3" = "green",
                               "C4" = "orange",
                               "C5" = "red"),
                   Group = c("Uninfected" = "grey",
                                  "Infected Recovered" = "black",
                                  "Infected PASC" = "magenta"),
                  `Fatigue...malaise` = c("NA" = "grey",
                                          "No" = "white",
                                          "Yes" = "black"),
                   "Pulmonary" = c("NA" = "grey",
                                    "No" = "white",
                                    "Yes" = "black"),
                    "Cardiovascular" = c("NA" = "grey",
                                         "No" = "white",
                                         "Yes" = "black"),
                    "Gastrointestinal" = c("NA" = "grey",
                                            "No" = "white",
                                            "Yes" = "black"),
                    "Musculoskeletal" = c("NA" = "grey",
                                          "No" = "white",
                                          "Yes" = "black"),
                    "Neurologic" = c("NA" = "grey",
                                     "No" = "white",
                                     "Yes" = "black"),
                    `Any.mild.symptom` = c("NA" = "grey",
                                           "No" = "white",
                                          "Yes" = "black"),
                    "Tinnitus" = c("NA" = "grey",
                                    "No" = "white",
                                    "Yes" = "black"))

# define custom color palatte
colorPalette <- c("purple", "black", "yellow")
colorPalette <- colorRampPalette(colors = colorPalette)(100)

# plot
outFile <- "Figure 4A.pdf"
pdf(file = outFile, width = 10, height = 10)
        
pheatmap(mat = subOlink[rownames(matannot_row), corder],
        #mat = subOlink[rownames(matannot_row), ],
        color = colorPalette,
        breaks = c(min(subOlink),
                   seq(from = -1*limit,
                       to = limit,
                       length.out = 99),
                   max(subOlink)),
        cellwidth = 2,
        cellheight = 4,
        cluster_cols = FALSE,
        cluster_rows = FALSE,
        annotation_col = matannot_col,
        annotation_row = matannot_row,
        show_rownames = TRUE,
        show_colnames = FALSE,
        annotation_colors = ann_colors,
        border_color = NA,
        fontsize_row = 4,
        fontsize = 7)
dev.off()



### Box and jitter plot of selected cytokines, chemokines expression across cluster
selproteins <- c("IFNG", "IL12B", "EBI3_IL27", "CXCL9", "CXCL10", "CXCL11",  "TNF", "IL6", "CCL7")

plotDF <- Olink_rulein[selproteins, , drop = F] %>%
          rownames_to_column(var = "protein") %>%
          gather(PTID_visit, NPX, -protein) %>%
          mutate(clusterID = pivotDF$Cluster[match(PTID_visit,
                            table = pivotDF$PTID_visit)],
                 Group = pivotDF$Group[match(PTID_visit,
                                    table = pivotDF$PTID_visit)])

plotDF$protein <- factor(plotDF$protein, levels = selproteins)
plotDF$Group <- factor(plotDF$Group, levels = names(gCols))

bPlot <- ggplot(data = plotDF,
                mapping = aes(x = clusterID, y = NPX)) +
          geom_boxplot(outlier.shape = NA,
                       size = 0.3) +
          geom_point(position = position_jitterdodge(),
                     aes(color = factor(Group)),
                         #shape = 21,
                     width = 0.1,
                     size = 1.5) +
          scale_color_manual(values = gCols) +
          facet_wrap(~protein, scales = "free", ncol = 6) +
          labs(x = NULL,
               y = "Normalized Protein Expression",
               color = "Group") +
          theme_bw() +
          theme(strip.text = element_text(size = 12),
                axis.text.x = element_text(color = "black", size = 9),
                axis.ticks.x = element_blank(),
                axis.text.y = element_text(size = 12, color = "black"),
                axis.line = element_line(color = "grey"),
                panel.background = element_blank(),
                panel.grid = element_blank())
pdf(file = "Figure 3B 3E.pdf", width = 35, height = 35, useDingbats = F)
print(bPlot)
dev.off()


#### Line plots of proteins over time for longitudinal samples
selproteins <- c("IFNG", "IL12B", "EBI3_IL27", "CXCL9", "CXCL10", "CXCL11",  "TNF", "IL6", "CCL7", "SAMD9L", "MNDA", "DDX58", "LAMP3")

## show inflammatory PASC (as C4+C5), Non-inflammatory PASC (C2+C3) and Recovered's as C1+C2+C3
pivotDF3 <- tableS1 %>%
            filter(!Group == "Uninfected") %>%
            mutate(iGroup = ifelse(Cluster %in% c(4, 5),
                                                  "Inflammatory Group",
                                                  "Non-inflammatory Group"),
                   iGroup2 = paste(Group, iGroup, sep = "_")) %>%
            filter(!iGroup2 == "Infected Recovered_Inflammatory Group") %>%
            mutate(iGroup2 = ifelse(iGroup2 == "Infected Recovered_Non-inflammatory Group",
                                    "Infected Recovered",
                                    iGroup2))

keeplong <- names(table(pivotDF3$PTID))[table(pivotDF3$PTID) > 1]
pivotDF3 <- pivotDF3 %>% filter(PTID %in% keeplong)

df1 <- pivotDF3 %>% select(PTID, iGroup2) %>% unique()
table(df1$iGroup2)

# keep samples selected from pivotDF3
olink_long <- fullOlink[selproteins, pivotDF3$PTID_visit, drop = F] %>%
              rownames_to_column(var = "protein") %>%
              gather(PTID_visit, NPX, -protein) %>%
              inner_join(., pivotDF3, by = "PTID_visit")
olink_long$"Days PSO" <- as.numeric(olink_long$"Days PSO")
olink_long$protein <- factor(olink_long$protein, levels = selproteins)

# remove >274 days since thats the last timepoint with atleast 3 PTIDs, everything after had 1 or 2
olink_long_v1 <- olink_long %>% filter(protein == "IFNG")
rowSums(table(olink_long_v1$"Days PSO", olink_long_v1$iGroup2))

olink_long <- olink_long %>% filter(`Days PSO` <=274)

# status color
cols1 <- c("Infected Recovered" = "black",
           "Infected PASC_Non-inflammatory Group" = "blue",
           "Infected PASC_Inflammatory Group" = "red")

# line plot with smoothing curve
lPlot <-  ggplot(data = olink_long,
                 mapping = aes(x = `Days PSO`,
                               y = NPX,
                               fill = factor(iGroup2),
                               color = factor(iGroup2))) +
              stat_smooth(method = "loess",
                          alpha = 0.4,
                          se = TRUE,
                          size = 0.5) +
           geom_vline(xintercept = 60, linetype = "dashed") +
              scale_fill_manual(values = cols1) +
              scale_color_manual(values = cols1) +
              facet_wrap(~protein, scales = "free", ncol = 7) +
              labs(x = "Days: symptom onset to blood draw",
                   y = "Normalized Protein Expression", color = NULL) +
              theme_bw() +
              theme(strip.text = element_text(size = 12),
                    axis.text.x = element_text(size = 10, color = "black"),
                    axis.text.y = element_text(size = 11, color = "black"),
                    axis.line = element_line(color = "grey"),
                    panel.background = element_blank(),
                    panel.grid = element_blank())
pdf(file = "Figure 3C 3F 3H.pdf", width = 12, height = 4)
print(lPlot)
dev.off()

# line plot showing all PTID points
lPlot <-  ggplot(data = olink_long,
                    mapping = aes(x = `Days PSO`,
                                  y = NPX,
                                  group = factor(PTID))) +
              #geom_point(aes(color = factor(iGroup2)), size = 1) +
              geom_line(aes(color = factor(iGroup2)), size = 0.2) +
              scale_color_manual(values = cols1) +
              geom_vline(xintercept = 60, linetype = "dashed") +
              facet_wrap(~protein, scales = "free", ncol = 6) +
              labs(x = "Days: symptom onset to blood draw",
                   y = "Normalized Protein Expression", color = NULL) +
              theme_bw() +
              theme(strip.text = element_text(size = 12),
                    axis.text.x = element_text(size = 10, color = "black"),
                    axis.text.y = element_text(size = 11, color = "black"),
                    axis.line = element_line(color = "grey"),
                    panel.background = element_blank(),
                    panel.grid = element_blank())
pdf(file = "Fig S9.pdf", width = 15, height = 5, useDingbats = F)
print(lPlot)
dev.off()



#### line plots of significant selected modules over time
selmodules <- c("biocarta tid pathway",
                "wp type ii interferon signaling ifng",
                "pid il27 pathway",
                "wp il18 signaling pathway",
                "pid nfkappab canonical pathway",
                "reactome regulation of ifna signaling",
                "pid il1 pathway",
                "reactome tnf signaling")

sigmods2 <- sigmods %>% filter(name %in% selmodules)
gsLS2 <- sigmods2 %>% .$genes %>% strsplit(",")
names(gsLS2) <- sigmods2$name

mat_exp <- fullOlink[, pivotDF3$PTID_visit, drop = F]
set.seed(20220122)
ssgsea_res <- gsva(data.matrix(mat_exp),
                   gsLS2,
                   method = "ssgsea",
                   abs.ranking = FALSE,
                   min.sz = 2,
                   max.sz = 2000,
                   ssgsea.norm = T)
ssgsea_res <- t(scale(t(ssgsea_res)))

ssgsea_long <- ssgsea_res %>%
               as.data.frame() %>%
               rownames_to_column(var = "module") %>%
               gather(PTID_visit, ssgsea_score, -module) %>%
               inner_join(., pivotDF3, by = "PTID_visit")
ssgsea_long$"Days PSO" <- as.numeric(ssgsea_long$"Days PSO")

# remove >274 days since thats the last timepoint with atleast 3 PTIDs, everything after had 1 or 2
ssgsea_long_v1 <- ssgsea_long %>% filter(module == "biocarta tid pathway")
rowSums(table(ssgsea_long_v1$"Days PSO", ssgsea_long_v1$iGroup2))

ssgsea_long <- ssgsea_long %>% filter(`Days PSO` <=274)

# line plot with smoothing curve
lPlot <-  ggplot(data = ssgsea_long,
                 mapping = aes(x = `Days PSO`,
                               y = ssgsea_score,
                               fill = factor(iGroup2),
                               color = factor(iGroup2))) +
              stat_smooth(method = "loess",
                          alpha = 0.4,
                          se = TRUE,
                          size = 0.5) +
           geom_vline(xintercept = 60, linetype = "dashed") +
              scale_fill_manual(values = cols1) +
              scale_color_manual(values = cols1) +
              facet_wrap(~module, scales = "free") +
              labs(x = "Days: symptom onset to blood draw",
                   y = "ssGSEA score", color = NULL) +
              theme_bw() +
              theme(strip.text = element_text(size = 12),
                    axis.text.x = element_text(size = 10, color = "black"),
                    axis.text.y = element_text(size = 11, color = "black"),
                    axis.line = element_line(color = "grey"),
                    panel.background = element_blank(),
                    panel.grid = element_blank())
pdf(file = "Figure 3D 3G 3I.pdf", width = 12, height = 8)
print(lPlot)
dev.off()


# line plot showing all PTID points
lPlot <-  ggplot(data = ssgsea_long,
                    mapping = aes(x = `Days PSO`,
                                  y = ssgsea_score,
                                  group = factor(PTID))) +
              geom_line(aes(color = factor(iGroup2)), size = 0.2) +
              scale_color_manual(values = cols1) +
              geom_vline(xintercept = 60, linetype = "dashed") +
              facet_wrap(~module, scales = "free") +
              labs(x = "Days: symptom onset to blood draw",
                   y = "ssGSEA score", color = NULL) +
              theme_bw() +
              theme(strip.text = element_text(size = 12),
                    axis.text.x = element_text(size = 10, color = "black"),
                    axis.text.y = element_text(size = 11, color = "black"),
                    axis.line = element_line(color = "grey"),
                    panel.background = element_blank(),
                    panel.grid = element_blank())
pdf(file = "Figure S10.pdf", width = 15, height = 10, useDingbats = F)
print(lPlot)
dev.off()




############################################################
############ External dataset validation ############
#### Jim heath proteomics data
############################################################
jhDF <- read_excel("JimHeath_proteomics.xlsx", sheet = 2) %>%
        as.data.frame() %>%
        mutate(subjectID = ifelse(`Blood Draw` %in% NA,
                                  `Patient Subject ID`,
                                  `Blood Draw`))
jh_meta <- read_excel("JimHeath_proteomics_metadata.xlsx", sheet = 2) %>%                as.data.frame()
jh_symptoms <- read_excel("JimHeath_proteomics_metadata.xlsx", sheet = 4) %>%                as.data.frame() %>%
               column_to_rownames(var = "Study subject ID")
jh_symptomsDF <- data.frame(nsymptoms = rowSums(jh_symptoms, na.rm = T)) %>%
                 rownames_to_column(var = "subjectID") %>%
                 mutate(symptomatic = ifelse(nsymptoms == 0, "INCOV no",
                                                             "INCOV yes"))

jh_annot <- jhDF %>%
            select(`Patient Subject ID`,
                    subjectID,
                   `Blood Draw`,
                   `Healthy or INCOV`,
                    age,
                    sex,
                    BMI) %>%
            mutate(timepoint = gsub(".+-(.+)", "\\1", `Blood Draw`),
                   WOS = jh_meta$"Who Ordinal Scale"[match(`Blood Draw`,
                                            table = jh_meta$"Blood draw")],
                   DaysPSO = jh_meta$"Observation days since onset of symptoms"[match(`Blood Draw`,
                                        table = jh_meta$"Blood draw")],
                   DaysPSO = round(DaysPSO),
                   symptomatic = jh_symptomsDF$symptomatic[match(`Patient Subject ID`,
                                table = jh_symptomsDF$subjectID)],
                   symptomatic = ifelse(`Healthy or INCOV` == "Healthy", "Healthy", symptomatic),
                   symptomatic = ifelse(symptomatic %in% NA, "INCOV not reported", symptomatic))


# filter on COV dayspSO > 60: take first timepoint at >60
# remove PTID with no symtpom recording
d60_ptids <- jh_annot %>%
             filter(DaysPSO >= 60 &
                   !symptomatic == "INCOV not reported") %>%
             group_by(`Patient Subject ID`) %>%
             top_n(-1, DaysPSO) %>%
             .$subjectID

healthyIDs <- jh_annot %>%
              filter(`Healthy or INCOV` == "Healthy") %>%
              .$subjectID

# final matrix
jhmat <- jhDF %>%
         filter(subjectID %in% c(healthyIDs, d60_ptids)) %>%
         select(-`Patient Subject ID`,
                -`Blood Draw`,
                -`Healthy or INCOV`,
                -age,
                -sex,
                -BMI) %>%
         column_to_rownames(var = "subjectID")
jhmat <- as.data.frame(t(jhmat))

# check overlap of AIFI full olink panel list with JH list
jhpanel <- unique(gsub("_.+", "", rownames(jhmat)))
ovlp1 <- intersect(rownames(fullOlink), jhpanel)


# Patient annotation dataframe
matannot_col <- jh_annot %>%
                filter(subjectID %in% colnames(jhmat)) %>%
                select(subjectID, `Healthy or INCOV`, age, sex, BMI, DaysPSO, symptomatic)

# add acute WOS score
acutewos <- jh_annot %>%
            filter(`Healthy or INCOV` == "INCOV" & timepoint == "T1") %>%
            select(subjectID, WOS) %>%
            mutate(PTID = gsub("-.+", "", subjectID))
matannot_col <- matannot_col %>%
                mutate(PTID = gsub("-.+", "", subjectID)) %>%
                mutate(aWOS = acutewos$WOS[match(PTID,
                                                table = acutewos$PTID)]) %>%
                select(-PTID) %>%
                column_to_rownames(var = "subjectID")

matannot_col <- matannot_col[colnames(jhmat), , drop = F]
table(rownames(matannot_col) == colnames(jhmat))


#### check overlap of clusters 4 and 5 biomarker proteins with JH cohort
clusterproteins <- clDF2 %>%
                   filter(adjp_val < 0.05) %>%
                   filter(clusterID %in% c("H4", "H5")) %>%
                  .$protein %>%
                   unique()
ovlp2 <- intersect(clusterproteins, jhpanel)

# subset JH data on cluster4 and cluster 5 defining overlapping proteins
pattern <- paste(ovlp2, collapse = "|")
selp <- c(rownames(jhmat)[grep(pattern, rownames(jhmat))],
          rownames(jhmat)[grep("IL27", rownames(jhmat))])
jhmat_ovlp <- jhmat[selp, , drop = F]
jhmat_ovlp <- t(scale(t(jhmat_ovlp)))


# find clusters using kmeans
library("factoextra")
set.seed(202203082)
kmat2 <- kmeans(t(jhmat_ovlp),
                centers = 5,
                iter.max = 100,
                nstart = 100,
                algorithm = "Hartigan-Wong")
kclusters2 <- fviz_cluster(kmat2,
                           data = t(jhmat_ovlp),
                           geom = "point",
                           shape = 16,
                           show.clust.cent = FALSE,
                           palette = c("lightblue",
                                       "pink",
                                       "lightgreen",
                                       "orange",
                                       "tan",
                                       "darkgrey",
                                       "black",
                                       "magenta"),
                           ggtheme = theme_classic()) +
              geom_point(aes(color = matannot_col$symptomatic)) +
              ggtitle("K-means clustering 163 proteins") +
              theme(axis.text.x = element_text(size = 10, color = "black"),
                    axis.text.y = element_text(size = 10, color = "black"))
pdf(file = "Figure 5A.pdf", width = 8, height = 5, useDingbats = FALSE)
print(kclusters2)
dev.off()


# get clusters info
clusterDF <- data.frame(clusterID = sort(kmat2$cluster)) %>%
             rownames_to_column(var = "subjectID")
matannot_col$clusterID2 <- clusterDF$clusterID[match(rownames(matannot_col),
                                        table = clusterDF$subjectID)]

table(matannot_col$"Healthy or INCOV", matannot_col$clusterID2)
    

# pie charts of group porportions per cluster
matannot_col2 <- matannot_col %>%
                 rownames_to_column()  %>%
                 mutate(clusterID3 = ifelse(clusterID2 == 5, "A",
                                     ifelse(clusterID2 == 3, "B",
                                     ifelse(clusterID2 == 1, "C",
                                     ifelse(clusterID2 == 4, "D",
                                     ifelse(clusterID2 == 2, "E", NA)))))) %>%
                 column_to_rownames(var = "rowname")

plotDF <- as.data.frame(table(matannot_col2$symptomatic,
                              matannot_col2$clusterID3)) %>%
          spread(Var1, Freq) %>%
          mutate(pHealthy = (Healthy/(Healthy+`INCOV no`+`INCOV yes`))*100,
                 pINCOVn = (`INCOV no`/(Healthy+`INCOV no`+`INCOV yes`))*100,
                 pINCOVy = (`INCOV yes`/(Healthy+`INCOV no`+`INCOV yes`))*100) %>%
          select(Var2, pHealthy, pINCOVn, pINCOVy) %>%
          gather(group, prop, -Var2)

pCols <- c("pHealthy" = "darkgrey",
           "pINCOVn" = "black",
           "pINCOVy" = "magenta")

plotpie <- ggplot(data = plotDF,
                   mapping =  aes(x = "Var2",
                                  y = prop,
                                  fill = group)) +
            scale_fill_manual(values = pCols) +
            geom_bar(stat = "identity", width = 1, color = "white") +
            coord_polar("y", start = 0) +
            facet_wrap(~Var2, nrow = 1) +
            theme_bw() +
            theme(axis.text.x = element_text(size = 12, color = "black"),
                  axis.text.y = element_text(size = 12, color = "black"),
                  axis.line = element_line(color = "grey"),
                  panel.background = element_blank(),
                  panel.grid = element_blank())
pdf(file = "Figure 5B.pdf", width = 10, height = 5)
print(plotpie)
dev.off()



#### find markers that are different between the two groups : inf VS non-inf only PASC
gatherDF <- jhmat_ovlp %>%
            as.data.frame() %>%
            rownames_to_column(var = "protein") %>%
            gather(subjectID, value, -protein)
gatherDF1 <- merge(gatherDF,
                   matannot_col2,
                   by.x = "subjectID",
                   by.y = "row.names") %>%
            mutate(groups_to_compare = ifelse(clusterID3 == "E",
                                              "inf",
                                              "noninf"),
                   group_annot = paste(groups_to_compare,
                                       `Healthy or INCOV`,
                                       sep = "_"))
gatherDF2 <- gatherDF1 %>% filter(`Healthy or INCOV` == "INCOV")

# calculate median protein expression per group
medDF <- gatherDF2 %>%
         group_by(protein, groups_to_compare) %>%
         summarize(m = median(value)) %>%
         as.data.frame() %>%
         spread(groups_to_compare, m) %>%
         mutate(delta = inf - noninf,
                group = ifelse(inf > noninf,
                               "Higher in INCOV E",
                               "Higher in INCOV B,C,D")) %>%
         arrange(group)

# wilcox tests cluster E vs B,C,C INCOV
wDF <- gatherDF2 %>%
       group_by(protein) %>%
       do(p = wilcox.test(value ~ groups_to_compare, data = .)$p.value) %>%
       mutate(p = unlist(p)) %>%
       as.data.frame() %>%
       mutate(adjp = p.adjust(p, method = "BH")) %>%
       filter(adjp < 0.05) %>%
       arrange(adjp) %>%
       inner_join(., medDF, by = "protein")
sort(wDF$protein)


# box plots of cluster E VS B,C,D incov proteins of select proteins
selproteins <- c("TNF_INF", "IL12B_INF", "CCL7_INF", "CXCL11_INF", "CXCL10_INF", "IFNG_INF")
gatherDF3 <- gatherDF1 %>%
             filter(!`Healthy or INCOV` == "Healthy") %>%
             filter(protein %in% selproteins)

gatherDF3$group_annot <- factor(gatherDF3$group_annot,
                                levels = names(gCols))

porder <- selproteins
gatherDF3$protein <- factor(gatherDF3$protein, levels = porder)

# box plots
gCols <- c("noninf_Healthy" = "darkgrey",
           "noninf_INCOV" = "tan",
           "inf_INCOV" = "purple")

bPlot <- ggplot(data = gatherDF3,
                mapping = aes(x = group_annot, y = value)) +
          geom_boxplot(outlier.shape = NA,
                       size = 0.3) +
          geom_jitter(aes(color = group_annot),
                     width = 0.2,
                     size = 2) +
          scale_color_manual(values = gCols) +
          facet_wrap(~protein, scales = "free_y") +
          labs(x = NULL, y = "NPX") +
          theme_bw() +
          theme(strip.text = element_text(size = 10),
                axis.text.x = element_text(color = "black", size = 5),
                axis.ticks.x = element_blank(),
                axis.text.y = element_text(size = 12, color = "black"),
                axis.line = element_line(color = "grey"),
                panel.background = element_blank(),
                panel.grid = element_blank())
pdf(file = "Figure 5C.pdf", width = 30, height = 40, useDingbats = F)
print(bPlot)
dev.off()


### plot including healthy controls for supplemental figure
porder <- wDF$protein[1:10]
gatherDF4 <- gatherDF1 %>% filter(protein %in% porder)
gatherDF4$group_annot <- factor(gatherDF4$group_annot,
                                levels = names(gCols))
gatherDF4$protein <- factor(gatherDF4$protein, levels = porder)

bPlot <- ggplot(data = gatherDF4,
                mapping = aes(x = group_annot, y = value)) +
          geom_boxplot(outlier.shape = NA,
                       size = 0.3) +
          geom_jitter(aes(color = group_annot),
                     width = 0.2,
                     size = 2) +
          scale_color_manual(values = gCols) +
          facet_wrap(~protein, scales = "free_y") +
          labs(x = NULL, y = "NPX") +
          theme_bw() +
          theme(strip.text = element_text(size = 10),
                axis.text.x = element_text(color = "black", size = 5),
                axis.ticks.x = element_blank(),
                axis.text.y = element_text(size = 12, color = "black"),
                axis.line = element_line(color = "grey"),
                panel.background = element_blank(),
                panel.grid = element_blank())
pdf(file = "Figure S12.pdf", width = 30, height = 40, useDingbats = F)
print(bPlot)
dev.off()



### check difference in acute WOS between clusters
# grouped aWOS
matannot_col_v2 <- matannot_col2 %>%
                   rownames_to_column(var = "subjectID") %>%
                   mutate(groups_to_compare = ifelse(clusterID3 == "E",
                                                     "inf",
                                                     "noninf"),
                          group_annot = paste(groups_to_compare,
                                              `Healthy or INCOV`,
                                              sep = "_"),
                          aWOS_grouped = ifelse(aWOS %in% c("1", "1 or 2"),
                                                "<=2",
                                         ifelse(aWOS %in% c("3", "4", "5"),
                                                "3 to 5",
                                         ifelse(aWOS %in% c("6", "7"),
                                                "6 and 7",
                                         NA)))) %>%
                    filter(`Healthy or INCOV` == "INCOV")

tabDF <- table(matannot_col_v2$group_annot, matannot_col_v2$aWOS_grouped)
fisher.test(tabDF)

# stacked bar plot of proportions from above table
cSums <- colSums(tabDF)
plotDF <- as.data.frame(t(tabDF) / colSums(tabDF)) %>%
          rename(WOS = Var1, INCOV_group = Var2)
          
plotbar <- ggplot(data = plotDF,
                   mapping =  aes(x = WOS,
                                  y = Freq,
                                  fill = INCOV_group,
                                  label = round(Freq, 1))) +
            geom_bar(position = "stack", stat = "identity") +
            geom_text(size = 6,
                      position = position_stack(vjust = 0.5),
                      color = "white") +
            scale_fill_manual(values = gCols) +
            labs(x = "Acute WHO ordinal scale",
                 y = "Propotion of patients",
                 fill = "INCOV group") +
            theme_bw() +
            theme(axis.text.x = element_text(size = 12, color = "black"),
                  axis.text.y = element_text(size = 12, color = "black"),
                  axis.line = element_line(color = "grey"),
                  panel.background = element_blank(),
                  panel.grid = element_blank())
pdf(file = "Figure 5D.pdf", width = 4, height = 5)
print(plotbar)
dev.off()


#######################################################
##### Figure 6: Diagnostic panel
#######################################################

## show bpx plots for diagnostic markers (INF vs NONINF PASC)
selproteins <- c("CCL7", "CD40LG", "S100A12")

plotDF <- Olink_rulein[selproteins, , drop = F] %>%
          rownames_to_column(var = "protein") %>%
          gather(PTID_visit, NPX, -protein) %>%
          mutate(clusterID = pivotDF$Cluster[match(PTID_visit,
                            table = pivotDF$PTID_visit)],
                 Group = pivotDF$Group[match(PTID_visit,
                                    table = pivotDF$PTID_visit)]) %>%
          filter(Group == "Infected PASC" & protein %in% selproteins) %>%
          mutate(iGroup = ifelse(clusterID %in% c("C4", "C5"),
                                                  "Inflammatory Group",
                                                  "Non-inflammatory Group"),
                   iGroup2 = paste(Group, iGroup, sep = "_"))
plotDF$protein <- factor(plotDF$protein, levels = selproteins)

# box and jitter plot
cols1 <- c("Infected PASC_Non-inflammatory Group" = "blue",
            "Infected PASC_Inflammatory Group" = "red")
plotDF$iGroup2 <- factor(plotDF$iGroup2, levels = names(cols1))

bPlot <- ggplot(data = plotDF,
                mapping = aes(x = clusterID, y = NPX)) +
          geom_boxplot(outlier.shape = NA,
                       size = 0.3) +
          geom_point(position = position_jitterdodge(),
                     aes(color = factor(iGroup2)),
                     width = 0.1,
                     size = 1.5) +
          scale_color_manual(values = cols1) +
          facet_wrap(~protein, scales = "free", ncol = 6) +
          labs(x = NULL,
               y = "Normalized Protein Expression",
               color = "Group") +
          theme_classic() +
          theme(strip.text = element_text(size = 12),
                axis.text.x = element_text(color = "black", size = 9),
                axis.ticks.x = element_blank(),
                axis.text.y = element_text(size = 12, color = "black"),
                panel.background = element_blank(),
                panel.grid = element_blank())
pdf(file = "Figure 6A.pdf", width = 9, height = 3, useDingbats = F)
print(bPlot)
dev.off()

write.xlsx(plotDF, "FigurePDFs_Tables_Manuscript/Source data/Source Figure 6A.xlsx", col.names = TRUE, row.names = TRUE, append = FALSE)



## show line plots for diagnostic markers over time
pivotDF3 <- tableS1 %>% filter(Group == "Infected PASC")
keeplong <- names(table(pivotDF3$PTID))[table(pivotDF3$PTID) > 1]
pivotDF3 <- pivotDF3 %>% filter(PTID %in% keeplong)

# keep samples selected from pivotDF3
olink_long <- fullOlink[selproteins, pivotDF3$PTID_visit, drop = F] %>%
              rownames_to_column(var = "protein") %>%
              gather(PTID_visit, NPX, -protein) %>%
              inner_join(., pivotDF3, by = "PTID_visit") %>%
              mutate(`Days PSO` = as.numeric(`Days PSO`))
olink_long$protein <- factor(olink_long$protein, levels = selproteins)

# remove >274 days since thats the last timepoint with atleast 3 PTIDs, everything after had 1 or 2
olink_long <- olink_long %>% filter(`Days PSO` <=274)

# status color
cols1 <- c("Non-inflammatory group" = "blue", "Inflammatory group" = "red")

# line plot with smoothing curve
lPlot <-  ggplot(data = olink_long,
                 mapping = aes(x = `Days PSO`,
                               y = NPX,
                               fill = factor(igroup),
                               color = factor(igroup))) +
              stat_smooth(method = "loess",
                          alpha = 0.4,
                          se = TRUE,
                          size = 0.5) +
           geom_vline(xintercept = 60, linetype = "dashed") +
              scale_fill_manual(values = cols1) +
              scale_color_manual(values = cols1) +
              facet_wrap(~protein, scales = "free", ncol = 7) +
              labs(x = "Days: symptom onset to blood draw",
                   y = "Normalized Protein Expression", color = NULL) +
              theme_bw() +
              theme(strip.text = element_text(size = 12),
                    axis.text.x = element_text(size = 10, color = "black"),
                    axis.text.y = element_text(size = 11, color = "black"),
                    axis.line = element_line(color = "grey"),
                    panel.background = element_blank(),
                    panel.grid = element_blank())
pdf(file = "Figure 6B.pdf", width = 15, height = 5)
print(lPlot)
dev.off()


### plot probability score of the 3 panel diagnostic markers
pscoreDF_aifi <- read.csv("AIFI_PASC_Olink_with_Score.csv") %>%
                 mutate(visit = gsub(".+_(.+)", "\\1", PTID_visit),
                        PUBID = pubids$PUBID[match(PTID, table = pubids$PTID)]) %>%
                 select(-PTID_visit) %>%
                 mutate(PTID_visit = paste(PUBID, visit, sep = "_"))
write.csv(pscoreDF_aifi, "code/AIFI_PASC_Olink_with_Score.csv")

pscoreDF_aifi <- read.csv("code/AIFI_PASC_Olink_with_Score.csv") %>%
                 select(PTID_visit, Group, iGroup, Score) %>%
                 mutate(Group_iGroup = paste(Group, iGroup, sep = "_")) %>%
                 filter(!Group_iGroup %in% c("Uninfected_Inflammatory Group",
                                            "Infected Recovered_Inflammatory Group"))
pscoreDF_pasc <- pscoreDF_aifi %>% filter(Group == "Infected PASC")
icols1 <- c("Non-inflammatory Group" = "blue", "Inflammatory Group" = "red")
icols2 <- c("Uninfected_Non-inflammatory Group" = "grey",
            "Infected Recovered_Non-inflammatory Group" = "black",
            "Infected PASC_Non-inflammatory Group" = "blue",
            "Infected PASC_Inflammatory Group" = "red")

pscoreDF_pasc$iGroup <- factor(pscoreDF_pasc$iGroup, levels = names(icols1))
bPlot <- ggplot(data = pscoreDF_pasc,
                mapping = aes(x = iGroup, y = Score)) +
         geom_boxplot(outlier.shape = NA, size = 0.3) +
          geom_jitter(aes(color = iGroup),
                     width = 0.2,
                     size = 2) +
          scale_color_manual(values = icols1) +
          labs(x = NULL, y = "Probability score") +
          theme_classic() +
          theme(strip.text = element_text(size = 10),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.text.y = element_text(size = 12, color = "black"),
                panel.background = element_blank(),
                panel.grid = element_blank())
pdf(file = "Figure 6D Training.pdf", width = 5, height = 3, useDingbats = F)
print(bPlot)
dev.off()


# plot with uninfected and recovered
pscoreDF_aifi$Group_iGroup <- factor(pscoreDF_aifi$Group_iGroup,
                                     levels = names(icols2))
bPlot <- ggplot(data = pscoreDF_aifi,
                mapping = aes(x = Group_iGroup, y = Score)) +
         geom_boxplot(outlier.shape = NA, size = 0.3) +
          geom_jitter(aes(color = Group_iGroup),
                     width = 0.2,
                     size = 2) +
          scale_color_manual(values = icols2) +
          labs(x = NULL, y = "Probability score") +
          theme_classic() +
          theme(strip.text = element_text(size = 10),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.text.y = element_text(size = 12, color = "black"),
                panel.background = element_blank(),
                panel.grid = element_blank())
pdf(file = "Figure S13 Training.pdf", width = 7, height = 3, useDingbats = F)
print(bPlot)
dev.off()


##### probability score on test data
pscoreDF_jh <- read.csv("JH_PASC_Olink_with_Score.csv") %>%
               select(Subject, `Healthy.or.INCOV`, symptomatic, clusterID3, Score, Status)
pscoreDF_pasc <- pscoreDF_jh %>%
                 filter(Status %in%
                        c("Inflammatory PASC", "Non-inflammatory PASC"))
icols1 <- c("Non-inflammatory PASC" = "blue", "Inflammatory PASC" = "red")
icols2 <- c("Healthy" = "grey",
            "No symptoms" = "black",
            "Non-inflammatory PASC" = "blue",
            "Inflammatory PASC" = "red")

pscoreDF_pasc$Status <- factor(pscoreDF_pasc$Status, levels = names(icols1))
bPlot <- ggplot(data = pscoreDF_pasc,
                mapping = aes(x = Status, y = Score)) +
          geom_boxplot(outlier.shape = NA, size = 0.3) +
          geom_jitter(aes(color = as.factor(Status)),
                      width = 0.3,
                      size = 2) +
          scale_color_manual(values = icols1) +
          labs(x = NULL, y = "Probability score") +
          theme_classic() +
          theme(strip.text = element_text(size = 10),
                axis.text.x = element_text(size = 8, color = "black"),
                axis.ticks.x = element_blank(),
                axis.text.y = element_text(size = 12, color = "black"),
                panel.background = element_blank(),
                panel.grid = element_blank())
pdf(file = "Figure 6D Test.pdf", width = 4, height = 3, useDingbats = F)
print(bPlot)
dev.off()



## plot with healthy and recovered
pscoreDF_jh$Status <- factor(pscoreDF_jh$Status, levels = names(icols2))
bPlot <- ggplot(data = pscoreDF_jh,
                mapping = aes(x = Status, y = Score)) +
          geom_boxplot(outlier.shape = NA, size = 0.3) +
          geom_jitter(aes(color = Status),
                     width = 0.2,
                     size = 2) +
          scale_color_manual(values = icols2) +
          labs(x = NULL, y = "Probability score") +
          theme_classic() +
          theme(strip.text = element_text(size = 10),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.text.y = element_text(size = 12, color = "black"),
                panel.background = element_blank(),
                panel.grid = element_blank())
pdf(file = "Figure S13 Test.pdf", width = 7, height = 3, useDingbats = F)
print(bPlot)
dev.off()



### calculate probability scores acorss all longitudinal samples
### performance on full data set >= 60 days PSO : PASC samples
pivotDF4 <- tableS1 %>%
            filter(Group == "Infected PASC") %>%
            mutate(`Days PSO` = as.numeric(`Days PSO`)) %>%
            mutate(daysbin2 = ifelse(`Days PSO` >=100 & `Days PSO` < 150,
                                       ">=100 & <150 Days PSO",
                                ifelse(`Days PSO` >=150 & `Days PSO` < 200,
                                       ">=150 & <200 Days PSO",
                                ifelse(`Days PSO` >=200,
                                      ">=200 Days PSO",
                                        `Days bin`)))) %>%
             mutate(daysbin2 = ifelse(daysbin2 == ">=60 Days PSO",
                                     ">=60 & < 100 Days PSO",
                                      daysbin2)) %>%
             filter(`Days PSO` >= 60) %>%
             mutate(igroup = ifelse(igroup == "Inflammatory group",
                                    "Inflammatory PASC",
                                    "Non-inflammatory PASC"))

# subset Olink data on the samples above and on the 3 protein diagnostic marker panel - refer to the "Diagnostic panel" R codes
proteins <- c("CCL7", "CD40LG", "S100A12")
tdata_long <- t(fullOlink[proteins, pivotDF4$PTID_visit, drop = F]) %>%
              as.data.frame() %>%
              rownames_to_column(var = "PTID_visit") %>%
              inner_join(., pivotDF4, by = "PTID_visit")

# performance on training
source("Diaonostic panel LogisticRegression.R")

bmodel <- c(-6.319462, 2.538293, 1.091321, 0.9560707) # refer to Supplementary table S10
names(bmodel) <- c("(Intercept)", "CCL7", "CD40LG", "S100A12")
tdata_long$Score <- addModelScore(tdata_long, model = bmodel)

tdata_long$daysbin2 <- factor(tdata_long$daysbin2,
                            levels = c(">=60 & < 100 Days PSO",
                                       ">=100 & <150 Days PSO",
                                        ">=150 & <200 Days PSO",
                                        ">=200 Days PSO"))
icols1 <- c("Non-inflammatory PASC" = "blue", "Inflammatory PASC" = "red")

tdata_long$igroup <- factor(tdata_long$igroup, levels = names(icols1))

bPlot <- ggplot(data = tdata_long,
                mapping = aes(x = daysbin2, y = Score)) +
         geom_boxplot(outlier.shape = NA, size = 0.3) +
         geom_point(aes(color = factor(igroup))) +
         geom_line(aes(group = PTID, color = factor(igroup)), size = 0.2,
                    linetype = "dashed") +
          scale_color_manual(values = icols1) +
          facet_wrap(~igroup) +
          labs(x = "Days PSO bin", y = "Probability score", color = "Group") +
          theme_classic() +
          theme(strip.text = element_text(size = 10),
                axis.text.x = element_text(size = 8, angle = 90, hjust = 1,
                                            color = "black"),
                axis.ticks.x = element_blank(),
                axis.text.y = element_text(size = 12, color = "black"),
                panel.background = element_blank(),
                panel.grid = element_blank())
pdf(file = "Figure S13B.pdf", width = 9, height = 6, useDingbats = F)
print(bPlot)
dev.off()


# wilcox tests
wdf <- tdata_long %>% filter(daysbin2 == ">=60 & < 100 Days PSO")
wilcox.test(Score~igroup, data = wdf)

