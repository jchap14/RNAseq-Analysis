##################### generate a median normalized matrix for heatmaps in R ############# 
R
# source("http://bioconductor.org/biocLite.R")
# biocLite("EBSeq")

#import data
b <- read.delim('Xinguo_RNAseq.raw_counts.matrix', quote='', check.names= FALSE)
a <- b
a[,1] = NULL
a <- as.matrix(a)
rownames(a) <- b[,1]
Sizes = EBSeq::MedianNorm(a)
# calculate 
NormData = EBSeq::GetNormalizedMat(a, Sizes)
#export normalized table
write.table(NormData, 'Xinguo_RNAseq.norm_counts.matrix', sep='\t', row.names= TRUE, col.names= TRUE, quote = FALSE)




