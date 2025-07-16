setwd
library(limma)
library(oncoPredict)
library(parallel)
library(limma)
library(ggplot2)
library(ggpubr)
library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)
set.seed(123)
data=read.table("FPKM.txt", header=T, sep="\t", check.names=F,row.names = 1,quote = "")
dimnames=list(rownames(data), colnames(data))
data=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames=dimnames)
colnames(data)=gsub("(.*?)\\_(.*?)", "\\2", colnames(data))
GDSC2_Expr=readRDS(file='GDSC2_Expr.rds')
GDSC2_Res=readRDS(file = 'GDSC2_Res.rds')
GDSC2_Res=exp(GDSC2_Res) 
calcPhenotype(trainingExprData = GDSC2_Expr,    
              trainingPtype = GDSC2_Res,        
              testExprData = data,              
              batchCorrect = 'eb',  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,    
              minNumSamples = 10,    
              printOutput = TRUE,    
              removeLowVaringGenesFrom = 'rawData')