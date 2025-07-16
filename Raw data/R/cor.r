setwd
library(Hmisc)
input_data <- read.csv("cor.csv", row.names = 1, check.names = F)
res<-rcorr(as.matrix(input_data))
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
a<-flattenCorrMatrix(res$r, res$P)
write.csv(a,"cor_result.csv")