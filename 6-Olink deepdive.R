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

# load sample annotation (Uninfected, first timepoint >=60 days PSO for PASC and last timepoint >60 PSO for recoverd) used for Rule-in method
pivotDF <- read_excel("c2-Clustering-FTP-Longitudinal.xlsx", sheet = 1) %>%
           as.data.frame() %>%
           mutate(iGroup = ifelse(Heatmap %in% c("H4", "H5"),
                                                    "Inflammatory Group",
                                                    "Non-inflammatory Group"))


## Clustering based on symptomology : only showing PASC
pivotDF2 <- read_excel("c2-Clustering-FTP-Longitudinal.xlsx", sheet = 2) %>%
            as.data.frame() %>%
            filter(Group == "Infected PASC" & daysbin == ">=60 Days PSO")
symptomsDF <- pivotDF2[, c(4,12:21)] %>%
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
matannot_col <- pivotDF2 %>%
                select(PTID, Sex) %>%
                unique()
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
p<- pheatmap(mat = symptomsDF,
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
nsamplecluster <- as.data.frame(table(pivotDF$Heatmap))
plotDF <- as.data.frame(table(pivotDF$Heatmap, pivotDF$Group)) %>%
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
pdf(file = "Fig S3A.pdf", width = 5, height = 4)
print(plotbar)
dev.off()



### Severity score difference between groups (ONLY AMONG PASC)
severityDF <- read_excel("Clinical activity score.xlsx", sheet = 1) %>%
              as.data.frame() %>%
              filter(PTID %in% pivotDF$PTID) %>%
              mutate(Group = pivotDF$Group[match(PTID,
                                                 table = pivotDF$PTID)],
                    clusterID = pivotDF$Heatmap[match(PTID,
                                                    table = pivotDF$PTID)],
                    iGroup = ifelse(clusterID %in% c("H4", "H5"),
                                                     "Inflammatory Group",
                                                     "Non-inflammatory Group")) %>%
              filter(Group == "Infected PASC")

# wilcox tests
wilcox.test(`Clinical activity score` ~ iGroup, data = severityDF)

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
pdf(file = "Fig 1B.pdf", width = 6, height = 6, useDingbats = F)
print(bPlot)
dev.off()


### RBD IgG estimates
abestDF <- read.csv("RBD IgG estimates.csv") %>%
         filter(`COVID.PTID` %in% pivotDF$PTID) %>%
         mutate(clusterID = pivotDF$Heatmap[match(`COVID.PTID`,
                           table = pivotDF$PTID)],
                Group = pivotDF$Group[match(`COVID.PTID`,
                               table = pivotDF$PTID)],
                iGroup = ifelse(clusterID %in% c("H4", "H5"),
                                    "Inflammatory Group",
                                    "Non-inflammatory Group")) %>%
        filter(!Group == "Uninfected")
# wilcox tests
wilcox.test(log10(IgGRBDd60_origscale) ~ iGroup, data = abestDF)

# plot boxplot
gCols <- c("Infected Recovered" = "black", "Infected PASC" = "magenta")
abestDF$Group <- factor(abestDF$Group, levels = names(gCols))
abestDF$clusterID <- gsub("H", "", abestDF$clusterID)

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
pdf(file = "Fig 1C.pdf", width = 10, height = 5, useDingbats = F)
print(bPlot)
dev.off()




## read olink data
load("Olink_NPX.Rda")

# keep samples selected from pivotDF
fullOlink <- mat2021[, pivotDF$PTID_visit, drop = F]

### ssGSEA on modules and rank modules by top expressed in each cluster
# read modules that define clustering
sigmods <- read_excel("c2-Clustering-FTP-Longitudinal.xlsx", sheet = 3) %>%
           as.data.frame() %>%
           mutate(name = tolower(name), name = gsub("_" ," ", name))
gsLS <- sigmods %>% .$genes %>% strsplit(",")
names(gsLS) <- sigmods$module

# run ssGSEA
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
clusterIDs <- unique(pivotDF$Heatmap)

clLS <- lapply(clusterIDs, function(n) {
           print(n)
                                        
           # extract samples that belong to the cluster
            p1samples <- pivotDF %>% filter(Heatmap == n) %>% .$PTID_visit
                                              
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
        
# write table
fileName <- "Cluster markers_modules.txt"
write.table(clDF1,
            file = fileName,
            row.names = FALSE,
            sep = "\t",
            quote = FALSE)
         


# network of signifincat modules that sigificantly define groups
graphDF <- clDF1 %>%
           filter(adjp_val < 0.05) %>%
           mutate(lp = -log10(adjp_val),
                  clusterID = gsub("H", "Cluster", clusterID)) %>%
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
outFile <- "Fig S4.pdf"
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

plotDF <- ssgsea_res[selModules, , drop = F] %>%
          as.data.frame() %>%
          rownames_to_column(var = "module") %>%
          gather(PTID_visit, ssgsea, -module) %>%
          mutate(clusterID = pivotDF$Heatmap[match(PTID_visit,
                            table = pivotDF$PTID_visit)],
                 clusterID = gsub("H", "C", clusterID),
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
                      size = 2) +
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
pdf(file = "Fig 1D_1E_1F_S5.pdf", width = 35, height = 23, useDingbats = F)
print(bPlot)
dev.off()




# Correlellogram on top modules defning C4 and C5
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
subs_sgsea_res <- subs_sgsea_res[selmodules, , drop = F]

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
pdf(file = "Fig 1G.pdf", width = 5, height = 4)
print(cplot)
dev.off()




######## Identify individual DEPs of each sub-cluster
clusterIDs <- unique(pivotDF$Heatmap)

clLS <- lapply(clusterIDs, function(n) {
           print(n)
                                        
           # extract samples that belong to the cluster
            p1samples <- pivotDF %>% filter(Heatmap == n) %>% .$PTID_visit
                                              
            # extract samples that belong all other clusters
            p2samples <- setdiff(pivotDF$PTID_visit, p1samples)
            
          # calculate delta of expression
          # difference in each protein expression between the selected cluster samples and all other samples
          olinkDF <- t(fullOlink)
                
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
        
# write table
fileName <- "Cluster markers DEPs.txt"
write.table(clDF2,
            file = fileName,
            row.names = FALSE,
            sep = "\t",
            quote = FALSE)
    



### Network diagram of cluster specific proteins
# graph data frame to plot cytokines and cytokine receptors
# select cytokines and cytokine receptors
selMarkers <- rownames(fullOlink)[grep("^IL|EBI3|^CCL|^CXC|IFN|^TNF|^C2|^C4BPB|CSF1|CSF3|TGFB1$|TGFBR2|TGFBR3", rownames(fullOlink))]
selMarkers <- selMarkers[-grep("RN|ILKAP|C2CD2L", selMarkers)]

# make graph data frame
graphDF <- clDF2 %>%
           filter(adjp_val < 0.05) %>%
           filter(protein %in% selMarkers) %>%
           mutate(lp = -log10(p_val),
                  clusterID = gsub("H", "Cluster", clusterID)) %>%
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
outFile <- "Fig 2A.pdf"
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
outFile <- "Fig S6.pdf"
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

plotDF <- fullOlink[selproteins, , drop = F] %>%
          rownames_to_column(var = "protein") %>%
          gather(PTID_visit, NPX, -protein) %>%
          mutate(clusterID = pivotDF$Heatmap[match(PTID_visit,
                            table = pivotDF$PTID_visit)],
                 Group = pivotDF$Group[match(PTID_visit,
                                    table = pivotDF$PTID_visit)])

plotDF$protein <- factor(plotDF$protein, levels = selproteins)

# box and jitter plot
gCols <- c("Uninfected" = "darkgrey",
          "Infected Recovered" = "black",
          "Infected PASC" = "magenta")
plotDF$Group <- factor(plotDF$Group, levels = names(gCols))
plotDF$clusterID <- as.factor(gsub("H", "", plotDF$clusterID))

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
pdf(file = "Fig 2B_2E_S7.pdf", width = 35, height = 35, useDingbats = F)
print(bPlot)
dev.off()




#### Line plots of proteins over time for longitudinal samples
selproteins <- c("IFNG", "IL12B", "EBI3_IL27", "CXCL9", "CXCL10", "CXCL11",  "TNF", "IL6", "CCL7", "SAMD9L", "MNDA", "DDX58", "LAMP3")

## show inflammatory PASC (from H5, H5), Non-inflammatory PASC (H2, H3) and Recovered's from C1,C2,C3

pivotDF3 <- read_excel("c2-Clustering-FTP-Longitudinal.xlsx", sheet = 2) %>%
            as.data.frame() %>%
            filter(!Group == "Uninfected") %>%
            select(-Heatmap) %>%
            mutate(clusterID = pivotDF$Heatmap[match(PTID,
                                table = pivotDF$PTID)],
                    iGroup = ifelse(clusterID %in% c("H4", "H5"),
                                                  "Inflammatory Group",
                                                  "Non-inflammatory Group"),
                   iGroup2 = paste(Group, iGroup, sep = "_"),
                   DaysPSO = daysDF$DaysPSO[match(PTID_visit,
                                        table = daysDF$PTID_visit)]) %>%
            #filter(DaysPSO >= 60) %>%
            filter(!iGroup2 == "Infected Recovered_Inflammatory Group") %>%
            mutate(iGroup2 = ifelse(iGroup2 == "Infected Recovered_Non-inflammatory Group",
                                    "Infected Recovered",
                                    iGroup2))

keeplong <- names(table(pivotDF3$PTID))[table(pivotDF3$PTID) > 1]
pivotDF3 <- pivotDF3 %>% filter(PTID %in% keeplong)

df1 <- pivotDF3 %>% select(PTID, iGroup2) %>% unique()
table(df1$iGroup2)

# keep samples selected from pivotDF3
olink_long <- mat2021[selproteins, pivotDF3$PTID_visit, drop = F] %>%
              rownames_to_column(var = "protein") %>%
              gather(PTID_visit, NPX, -protein) %>%
              inner_join(., pivotDF3, by = "PTID_visit")

olink_long$protein <- factor(olink_long$protein, levels = selproteins)

# remove >274 days since thats the last timepoint with atleast 3 PTIDs, everything after had 1 or 2
olink_long_v1 <- olink_long %>% filter(protein == "IFNG")
rowSums(table(olink_long_v1$DaysPSO, olink_long_v1$iGroup2))

olink_long <- olink_long %>% filter(DaysPSO <=274)

# status color
cols1 <- c("Infected Recovered" = "black",
            "Infected PASC_Non-inflammatory Group" = "blue",
            "Infected PASC_Inflammatory Group" = "red")

# line plot with smoothing curve
lPlot <-  ggplot(data = olink_long,
                 mapping = aes(x = DaysPSO,
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
pdf(file = "Fig 2C_2F_2H.pdf", width = 15, height = 5)
print(lPlot)
dev.off()


# line plot showing all PTID points
lPlot <-  ggplot(data = olink_long,
                    mapping = aes(x = DaysPSO,
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
pdf(file = "Fig S8.pdf", width = 15, height = 5, useDingbats = F)
print(lPlot)
dev.off()



#### line plots of significant modules over time
sigmods2 <- sigmods %>% filter(name %in% selmodules)
gsLS2 <- sigmods2 %>% .$genes %>% strsplit(",")
names(gsLS2) <- sigmods2$name

mat_exp <- mat2021[, pivotDF3$PTID_visit, drop = F]
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

# remove >274 days since thats the last timepoint with atleast 3 PTIDs, everything after had 1 or 2
ssgsea_long_v1 <- ssgsea_long %>% filter(module == "biocarta tid pathway")
rowSums(table(ssgsea_long_v1$DaysPSO, ssgsea_long_v1$iGroup2))

ssgsea_long <- ssgsea_long %>% filter(DaysPSO <=274)

# line plot with smoothing curve
lPlot <-  ggplot(data = ssgsea_long,
                 mapping = aes(x = DaysPSO,
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
pdf(file = "Fig 2D_2G_2I.pdf", width = 12, height = 8)
print(lPlot)
dev.off()


# line plot showing all PTID points
lPlot <-  ggplot(data = ssgsea_long,
                    mapping = aes(x = DaysPSO,
                                  y = ssgsea_score,
                                  group = factor(PTID))) +
              #geom_point(aes(color = factor(iGroup2)), size = 1) +
              geom_line(aes(color = factor(iGroup2)), size = 0.2) +
              scale_color_manual(values = cols1) +
              #scale_x_discrete(breaks = breaks) +
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
pdf(file = "Fig S9.pdf", width = 15, height = 10, useDingbats = F)
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




#### check overlap of clusters 4 and 5 biomarker proteins
clusterproteins <- clDF2 %>%
                   filter(adjp_val < 0.05) %>%
                   filter(clusterID %in% c("H4", "H5")) %>%
                   .$protein %>%
                    unique()
ovlp2 <- c(intersect(clusterproteins, jhpanel), "EBI3_IL27")


# subset JH data on cluster4 and cluster 5 defining overlapping proteins
pattern <- paste(ovlp2, collapse = "|")
selp <- c(rownames(jhmat2)[grep(pattern, rownames(jhmat))],
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
pdf(file = "Fig 2J.pdf", width = 8, height = 5, useDingbats = FALSE)
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
pdf(file = "Fig 2J piechart.pdf", width = 10, height = 5)
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


# box plots of cluster E VS B,C,D incov proteins
gatherDF3 <- gatherDF1 %>%
             #filter(!`Healthy or INCOV` == "Healthy") %>%
             filter(protein %in% wDF$protein)

# box plots
gCols <- c("noninf_Healthy" = "darkgrey",
           "noninf_INCOV" = "tan",
           "inf_INCOV" = "purple")

gatherDF3$group_annot <- factor(gatherDF3$group_annot,
                                levels = names(gCols))

porder <- wDF$protein
gatherDF3$protein <- factor(gatherDF3$protein, levels = porder)

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
pdf(file = "Fig 2K_S10.pdf", width = 30, height = 40, useDingbats = F)
print(bPlot)
dev.off()


# check difference in acute WOS between clusters
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
pdf(file = "Fig 2L.pdf", width = 4, height = 5)
print(plotbar)
dev.off()
