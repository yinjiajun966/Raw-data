setwd
library(dplyr)
library(data.table)
library(GSVA)
library(Biobase)
library(stringr)
library(RColorBrewer)
library(pheatmap)
library(tidyr)  
library(ggpubr)
library(cowplot)
library(psych) 

ann=read.csv("anno.csv",header=TRUE, row.names=1,check.names = FALSE)
ann <- arrange(ann,desc(group))
gene_set<- read.csv('./EMT.csv',header = T)
gene_set<-gene_set[, 1:2]  
list<- split(as.matrix(gene_set)[,1], gene_set[,2]) 

exprSet <- read.csv("./FPKM.csv", header=TRUE,row.names=1,check.names = FALSE)
exprSet <- exprSet[,row.names(ann)]

gsva_matrix <- gsva(as.matrix(exprSet), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
write.table(gsva_matrix,file = 'ssGSEA_EMT.csv',quote = F,sep = ',')

