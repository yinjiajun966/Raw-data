setwd
library(e1071)
library(preprocessCore)
library(parallel)
library(CIBERSORT)
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)

data(LM22)
df<-read.csv("FPKM.csv",  row.names=1, header=T, check.names=F)
df=as.matrix(df)
results <- cibersort(sig_matrix = LM22, mixture_file = df,perm = 1000, QN = T)
write.csv(results,"cibersort.csv")

