setwd
library(WGCNA)
library(reshape2)
library(stringr)
library(data.table)
library(dplyr)
options(stringsAsFactors = FALSE)
enableWGCNAThreads(20)
dataExpr <- read.csv("FPKM.csv",  row.names=1, header=T, check.names=F)
dataExpr[,] = as.numeric(unlist(dataExpr[,]))  
dataExpr = dataExpr[complete.cases(dataExpr),]
dim(dataExpr)
head(dataExpr)[,1:8]
dataExpr <- as.data.frame(t(dataExpr))
gsg = goodSamplesGenes(dataExpr, verbose = 3)
gsg$allOK
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(dataExpr)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ", ")))
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(dataExpr), method = "complete")

pdf(file = "Fig1.Sample_cluster.pdf", width = 15, height = 5)
par(mar = c(0,4,2,0),xpd=F)
plot(sampleTree, main = "Sample clustering to detect outliers", 
     cex=0.8,sub="", 
     xlab="",
     cex.main=1.5)
dev.off()

dataTraits <- read.csv('trait.csv',header = T,check.names = F,row.names = 1)
enableWGCNAThreads(20)  
powers =seq(from = 1, to=10, by=1) 
sft = pickSoftThreshold(dataExpr, powerVector = powers,  verbose = 1) 
save(sft,file='sft.Rdata')

pdf('Fig2.Soft.Threshold.pdf',width = 10,height = 6)
par(mfrow = c(1,2))
cex1 = 0.9

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.90,col="red") 

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

softPower =3
adjacency = adjacency(dataExpr, power = softPower) 
TOM = TOMsimilarity(adjacency) 
dissTOM = 1-TOM 

geneTree = hclust(as.dist(dissTOM), method = "average")
pdf(file="Fig3.Gene.cluster.pdf",width=6,height=5)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

minModuleSize =  200 
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
pdf(file="Fig4.DynamicTree.pdf",width=6.5,height=4.5)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
MEList = moduleEigengenes(dataExpr, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average")

pdf('Fig5.merge.cluster.pdf',width = 9,height = 5)
par(mar=c(2,3,3,3))
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.2 
abline(h=MEDissThres, col = "red")
dev.off()

merge = mergeCloseModules(dataExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
pdf(file="Fig6.DynamicTree.pdf",width=6.5,height=5)
plotDendroAndColors(geneTree, cbind(dynamicColors,mergedColors),c("Dynamic Tree Cut","Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

moduleColors = mergedColors
save(moduleColors,file = 'moduleColors.Rdata')
table(moduleColors)
colorOrder = c("grey", standardColors(30))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
moduleTraitCor = cor(MEs, dataTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
write.table(moduleTraitPvalue,'moduleTraitPvalue.txt',sep = '\t',quote = F,row.names = T)
write.table(moduleTraitCor,'moduleTraitCor.txt',sep = '\t',quote = F,row.names = T)
pdf(file="Fig7.Module_trait.pdf",width=10,height=10)
par(mar = c(10, 10, 3, 10))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(dataTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(500),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text =0.8,
               main = paste("Module-trait relationships"))
dev.off()

modNames <- unique(moduleColors)
geneModuleMembership = as.data.frame(cor(dataExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
traitNames=colnames(dataTraits)
geneTraitSignificance = as.data.frame(cor(dataExpr, dataTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")

pdf('Fig8.Gene.trait.green.pdf',width = 8,height = 8)
par(mar=c(6,5,5,3))
for (i in colnames(geneTraitSignificance)){
  moduleGenes = moduleColors=='green'
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, 'MEgreen']),
                     abs(geneTraitSignificance[moduleGenes,i]),
                     xlab = ("Module Membership in green module"),
                     ylab = paste("Gene significance"),
                     main = paste0('Trait ',substr(i,4,nchar(i)),"\n"),
                     cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, col = 'green')
  abline(v=0.8,h=0.4,col="red")
}
dev.off()

probes = colnames(dataExpr)
geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)
for (Tra in 1:ncol(geneTraitSignificance))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[,Tra],
                         GSPvalue[, Tra])
  names(geneInfo0) = c(oldNames,names(geneTraitSignificance)[Tra],
                       names(GSPvalue)[Tra])
}

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,mod],
                         MMPvalue[, mod])
  names(geneInfo0) = c(oldNames,names(geneModuleMembership)[mod],
                       names(MMPvalue)[mod])
}
geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]
write.table(geneInfo, file = "GS_MM.csv",sep=",",row.names=F)

color_gene<-c("probes","moduleColor")

modulecolor_gene <- geneInfo[,colnames(geneInfo) %in% color_gene]
write.csv(modulecolor_gene,"modulegene.csv")
