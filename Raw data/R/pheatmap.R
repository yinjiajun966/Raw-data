setwd
library(pheatmap)          
rt<-read.csv("pheatmap.csv",header=T,row.names=1,check.names=F)    
group=read.csv("group.csv",header=T,row.names=1,check.names=F)   
rt<-log2(rt+1)
rt<-scale(rt,center=TRUE, scale=FALSE)
outFile="heatmap.pdf"      

pdf(file=outFile,width=6,height=5.5)
pheatmap(rt,
         annotation=group,
         cluster_cols = F,
         color = colorRampPalette(c("blue","skyblue","white","orange","red"))(50),
         show_colnames = F,
         scale="row",  
         border_color ="NA",
         fontsize = 8,
         fontsize_row=6,
         fontsize_col=6)
dev.off()