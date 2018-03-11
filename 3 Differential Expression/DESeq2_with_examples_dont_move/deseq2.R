#################################################################################################
######################### Differential Expression testing with DEseq2 ###########################
#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2"); biocLite("pasilla"); biocLite("DESeq")
library("DESeq2"); library("pasilla"); library("Biobase"); data("pasillaGenes")

##### Read in the count matrix
y <- read.delim('raw_counts.matrix.txt', quote='')
y$X <- make.names(y$X, unique = TRUE); rownames(y) <- y$X; y$X = NULL
y.matrix <- as.matrix(y); class(y.matrix) <- "integer"

##### Read in metadata file
metadata <- read.delim('metadata.txt', quote='')
rownames(metadata) <- metadata$X; metadata$X = NULL

#################################################################################################
##################### construct the DESeqDataSet with NO BLOCKING FACTORS #######################
dds <- DESeqDataSetFromMatrix(countData = y.matrix, colData = metadata, design = ~ Condition)
dds <- dds[ rowSums(counts(dds)) > 1, ] ## pre-filter to remove rows that have only 0 or 1 read
dds$condition <- relevel(dds$Condition, ref="Control") ## specify which factor is the control
##### Collapse technical replicates if you want (look up collapseReplicates function)
##### Differential Expression testing
dds <- DESeq(dds)
res <- results(dds)
summary(res) ## print summary of results
resOrdered <- res[order(res$padj),] ## order the results table by the smallest adjusted p value
resSig <- subset(resOrdered, padj < 0.1) # if wanting to subset results to a min padj
##### Exporting results
write.csv(as.data.frame(resOrdered), file="DEG_results_full.csv") 
write.csv(as.data.frame(resSig), file="DEG_results_FDR10.csv")
##### MA Plot
png('DEG_results.MAplot.png', width = 2000, height = 2000, res = 300, pointsize = 10)
plotMA(res, main="DESeq2", ylim=c(-2,2))
dev.off()
##### Plot counts for a single gene across the groups
plotCounts(dds, gene="ENSG00000215030.4", intgroup="condition")

#################################################################################################
##################### construct the DESeqDataSet with BLOCKING FACTORS #######################
ddsMF <- dds # create a copy of the DESeqDataSet, so that we can rerun the analysis using a multi-factor design.
design(ddsMF) <- formula(~ Gender + Condition) # first factor is the blocker
ddsMF <- DESeq(ddsMF) # rerun DEseq 
resMF <- results(ddsMF)
summary(resMF) ## print summary of results
resOrderedMF <- resMF[order(resMF$padj),] ## order the results table by the smallest adjusted p value
resSigMF <- subset(resOrderedMF, padj < 0.1) # if wanting to subset results to a min padj
##### Exporting results
write.csv(as.data.frame(resOrderedMF), file="DEG_resultsMF_full.csv") 
write.csv(as.data.frame(resSigMF), file="DEG_resultsMF_FDR10.csv")
##### MA Plot
png('DEG_resultsMF.MAplot.png', width = 2000, height = 2000, res = 300, pointsize = 10)
plotMA(resMF, main="DESeq2", ylim=c(-2,2))
dev.off()
##### Plot counts for a single gene across the groups
plotCounts(ddsMF, gene="ENSG00000215030.4", intgroup="condition")
