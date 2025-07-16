setwd
library('DESeq2')
countData <- read.csv("Count.csv",header=T,row.names=1,comment.char="",check.names=F)
colData <- read.csv("group.csv",header=T,comment.char="",check.names=F)
dds <- DESeqDataSetFromMatrix(countData = countData,colData = colData,design = ~ condition)
dds$condition <- relevel(dds$condition, ref = "control")
dds <- DESeq(dds)
res <- results(dds)
write.csv(res,"DESeq2.diff.csv")
