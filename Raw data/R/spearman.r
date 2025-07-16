setwd
library(psych)
library(stringr)
library(pheatmap)
df<-read.csv("cibersort_cor.csv",header=T,row.names = 1,check.names = F,stringsAsFactors = F)
df2<-read.csv("biomarker.csv",header=T,row.names = 1,check.names = F,stringsAsFactors = F)
a<-corr.test(df,df2,method = "spearman",adjust = "fdr")
cmt<-a$r
pmt<-a$p
if (!is.null(pmt)){
  sssmt<- pmt< 0.001
  pmt[sssmt] <-'***'
  ssmt <- pmt< 0.01 & pmt >0.001
  pmt[ssmt] <-'**'
  smt <- pmt >0.01& pmt <0.05
  pmt[smt] <- '*'
  pmt[!sssmt&!ssmt&!smt]<- ''
} else {
  pmt <- F
}
a<-pheatmap(cmt,scale = "none",cluster_row = F, cluster_col = F, border=NA,
            display_numbers = pmt,fontsize_number = 12, number_color = "black",
            cellwidth = 20, cellheight =20)

pdf('corheatmap_cibersort.pdf',width = 8,height = 8)
a
dev.off()