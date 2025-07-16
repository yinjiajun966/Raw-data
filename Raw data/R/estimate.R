setwd
library(estimate)
filterCommonGenes(input.f="FPKM.txt", output.f="estimate.gct", id="GeneSymbol")
estimateScore(input.ds="estimate.gct", output.ds="estimate_score.gct", platform="agilent")
estimate_score <- read.table("estimate_score.gct", skip = 2, header = TRUE)
write.csv(estimate_score,"estimate_est.csv",row.names = FALSE)