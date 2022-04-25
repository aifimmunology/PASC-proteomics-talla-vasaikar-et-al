# RuleIn criteria based identified pathways clubbed into modules
# convert these modules in GMT format

# cleanup & set directors
rm(list=ls())

################################################################################
# Calculate Jaccard index of identified modules
# Regenerate using Modules identified
c2_cluster <- read.csv("C2-Biomarker-pathway-Cluster.csv", stringsAsFactors = F)
c2_cluster <- c2_cluster[c2_cluster$Biomarker %in% "Yes",]$Marker

library(igraph)
gmxFile <- "c2.cp.v7.2.symbols.gmt"
colNames <- max(count.fields(file = gmxFile, sep = "\t"))
colNames <- seq(from = 1, to = colNames)
colNames <- as.character(colNames)
gmx <- read.table(file = gmxFile,
                  sep = "\t",
                  quote = "\"",
                  fill = TRUE,
                  col.names = colNames,
                  row.names = 1)
gmx <- gmx[, -1]
#intersect biomarker pathways
gmx <- gmx[c2_cluster,]
gmx <- apply(gmx, MARGIN = 1, FUN = function(x) {
  return(value = setdiff(unname(x), ""))
})
names(gmx) <- toupper(names(gmx))

# calculate Jaccard index
gsDist <- sapply(gmx, function(gs1) {
  gsDist <- sapply(gmx, function(gs2) {
    gs1 <- gs1
    gs2 <- gs2
    jaccardIndex <- length(intersect(gs1, gs2))
    jaccardIndex <- jaccardIndex /
      (length(gs1) + length(gs2) - jaccardIndex)
    return(value = jaccardIndex)
  })
  return(value = gsDist)
})

# define jaccardCutoff (0.25 or 0.5)
jaccardCutoff = 0.25
gsMin <- gsDist
gsMin[gsMin < jaccardCutoff] <- 0

# create graph
g <- graph.adjacency(gsMin,
                     mode = "undirected",
                     weighted = TRUE,
                     diag = FALSE)

# decompose graph into modules and
gD <- decompose(g)

# group modules
groupLS <- sapply(gD, FUN = function(x) match(V(x)$name, table = V(g)$name))

## find genes of each module
outNet <- lapply(1:length(gD), FUN = function(i) {
  gsName <- V(gD[[i]])$"name"
  geneFreq <- sort(table(unlist(gmx[gsName])))
  topGS <- gsName
  geneFreq <- geneFreq[order(geneFreq, decreasing = TRUE)]
  return(c(gs = paste(topGS, collapse = ","),
           genes = paste(names(geneFreq), collapse = ","),
           counts = length(geneFreq),
           freq  = paste(geneFreq, collapse = ","),
           module = paste("module", i, sep = "")))
})
outNet <- do.call(what = rbind, outNet)
outNet <- outNet %>% data.frame() %>% mutate(counts=as.numeric(counts))
# write results
write.table(outNet, file = "xC2-PathwayModules-v1.txt", sep = "\t", quote = F, row.names = F)
hist(outNet$counts, breaks = 100)

### Network of modules
vCol = "red"

outFile <- "PathwayModules.pdf"
pdf(file = outFile, width = 8, height = 8)
plot(g,
     mark.groups = groupLS,
     #mark.col = "transparent",
     edge.width = E(g)$"weight" * 5,
     vertex.color = vCol,
     vertex.frame.color = NA,
     vertex.label.color = "black",
     vertex.label.cex = 0.4,
     vertex.size = 6)
legend(x = 0.5,
       y = -1.1,
       lwd = range(E(g)$"weight") * 5,
       col = "grey",
       legend = signif(range(E(g)$"weight"), digits = 2),
       bty = "n",
       xjust = 0.5,
       title = "Jaccard index")
dev.off()
rm(g, gD, gmx, groupLS, gsDist, gsMin, colNames, jaccardCutoff, gmxFile, outFile, vCol)
################################################################################
