setwd
library(rjson)
metadata <- fromJSON(file = "metadata.cart.2024-07-31.json")
nSam <- length(metadata)
file_name <- vapply(seq_len(nSam), function(i) gsub(metadata[[i]]$file_name,pattern = '\\.gz',replacement = ''), 
                    character(1))
file_id <- vapply(seq_len(nSam), function(i) metadata[[i]]$file_id, 
                  character(1))
submitter_id <- vapply(seq_len(nSam), function(i) metadata[[i]]$associated_entities[[1]]$entity_submitter_id, 
                       character(1))
sample <- vapply(submitter_id, function(v) substr(v, 1, nchar(v) - 
                                                    16), character(1))
entity_submitter_id <- vapply(seq_len(nSam), function(i) metadata[[i]]$associated_entities[[1]]$entity_submitter_id, 
                              character(1))
sample_type <- vapply(submitter_id, function(i) substr(i, 14, nchar(i) - 
                                                         13), character(1))
sample_type[which(sample_type == '01')] <- 'PrimaryTumor'
sample_type[which(sample_type == '11')] <- 'SolidTissueNormal'
sample_type <- as.factor(sample_type)
patient <- vapply(submitter_id, function(v) substr(v, 1, nchar(v) - 
                                                     12), character(1))
metaMatrix <- data.frame(file_name, file_id, patient, sample, 
                         submitter_id, entity_submitter_id, sample_type, stringsAsFactors = FALSE)
metaMatrix <- metaMatrix[order(metaMatrix$submitter_id),    ]
metaMatrix[metaMatrix == "not reported"] <- NA
metaMatrix$sample_type <- gsub(" ", "", metaMatrix$sample_type, 
                               fixed = TRUE)
metaMatrix <- metaMatrix[grep('augmented_star_gene_counts',metaMatrix$file_name),]



cleanMirFun <- function(fl) {
  expr <- read.table(fl, header=T, stringsAsFactors = FALSE,sep = '\t')
  expr <- expr[-(1:4),]
  genenames <- expr[,2]
  expr <- expr$unstranded
  names(expr) <- genenames
  return(expr)
}

mRNAMatrix <- lapply(paste('./gdc_download_20240731_070933.731613.tar',metaMatrix$file_id,metaMatrix$file_name,sep = '/'), function(fl) cleanMirFun(fl))

STAD <- do.call("cbind", mRNAMatrix)
colnames(STAD) <- metaMatrix$patient
table(metaMatrix$sample_type)

group <- ifelse(as.character(sapply(colnames(STAD),function(x) substr(x,14,15)))== '01','Tumor','Normal')
write.table(STAD,'TCGA-STAD expdata Count.txt',sep = '\t')


