setwd
library(UpSetR)
outFile="upset.txt"        
outPic="upset.pdf"                  

files=dir()                         
files=grep("txt$",files,value=T)    
geneList=list()

for(i in 1:length(files)){
  inputFile=files[i]
  if(inputFile==outFile){next}
  rt=read.table(inputFile,header=F)       
  geneNames=as.vector(rt[,1])              
  geneNames=gsub("^ | $","",geneNames)     
  uniqGene=unique(geneNames)                
  header=unlist(strsplit(inputFile,"\\.|\\-"))
  geneList[[header[1]]]=uniqGene
  uniqLength=length(uniqGene)
  print(paste(header[1],uniqLength,sep=" "))
}


upsetData=fromList(geneList)
pdf(file=outPic,onefile = FALSE,width=9,height=6)
upset(upsetData,
      nsets = length(geneList),              
      nintersects = 50,                     
      order.by = "freq",                   
      show.numbers = "yes",               
      number.angles = 20,                  
      point.size = 2,                        
      matrix.color="red",                     
      line.size = 0.8,                       
      mainbar.y.label = "Gene Intersections", 
      sets.x.label = "Set Size")
dev.off()

intersectGenes=Reduce(intersect,geneList)
write.table(file=outFile,intersectGenes,sep="\t",quote=F,col.names=F,row.names=F)