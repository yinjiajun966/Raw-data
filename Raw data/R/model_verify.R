setwd
library(survminer)
module_gene<-read.delim('gene.txt',sep='\t',header = F)
module_gene=as.character(module_gene$V1)
module_gene
GSE31210_exp<-read.delim('GSE15459_series_matrix.txt',sep='\t',header = T)
GSE31210_exp[1:4,1:4]

gene_anno<-read.delim('GPL570-55999.txt',sep='\t',header = T,check.names = F)
gene_anno=gene_anno[,c("ID","Gene Symbol")]
gene_anno[,3]=stringr::str_split_fixed(gene_anno$`Gene Symbol`, ' /// ', 3)[, 2]
gene_anno=gene_anno[which(gene_anno$V3==''),c(1,2)]
gene_anno=gene_anno[which(gene_anno$`Gene Symbol`!=''),]
head(gene_anno)
colnames(gene_anno)=c('ID_REF','symbol')
GSE31210_exp=merge(gene_anno,GSE31210_exp,by='ID_REF')
GSE31210_exp[1:4,1:4]
GSE31210_exp=GSE31210_exp[,-1]
GSE31210_exp=aggregate(.~symbol,GSE31210_exp,mean)
GSE31210_exp[1:4,1:4]
rownames(GSE31210_exp)=GSE31210_exp$symbol
GSE31210_exp=GSE31210_exp[,-1]
GSE31210_exp[1:4,1:4]

GSE31210_cli<-read.delim('clinical.txt',sep='\t',header = T,check.names = F)
head(GSE31210_cli)
colnames(GSE31210_cli)=c('Samples','OS','OS.time')
GSE31210_cli=na.omit(GSE31210_cli)
GSE31210_cli$OS=gsub('dead',1,GSE31210_cli$OS)
GSE31210_cli$OS=gsub('alive',0,GSE31210_cli$OS)

com_gene1=intersect(rownames(GSE31210_exp),module_gene)
com_gene1
GSE31210_module<-data.frame(time=GSE31210_cli$OS.time,
                            status=GSE31210_cli$OS,
                            t(GSE31210_exp[com_gene1,GSE31210_cli$Samples]))
head(GSE31210_module)
GSE31210_module$time=GSE31210_module$time/12
GSE31210_module$status=as.integer(GSE31210_module$status)

library(survival)
fmla <- as.formula(paste0("Surv(time, status) ~"
                          ,paste0(com_gene1,collapse = '+')))

cox1 <- coxph(fmla, data =as.data.frame(GSE31210_module))
lan <- coef(cox1)
round(lan, 5)
genes <- names(cox1$coefficients)
paste0(round(lan, 5), '*', names(lan))
risk.GSE31210 <- as.numeric(lan%*%as.matrix(t(GSE31210_module[,genes])))
risk.GSE31210z <- mosaic::zscore(risk.GSE31210)

library(pROC)
ROC_rt=timeROC::timeROC(T=GSE31210_module$time, 
                        delta=GSE31210_module$status,
                        marker=risk.GSE31210z, cause=1,
                        weighting='marginal',
                        times=c(1,3,5), 
                        ROC=TRUE,iid = T)
ROC_rt$AUC

pdf('geo_ROC.pdf',he=7,wi=7)
plot(ROC_rt, time=1, col='red', title=FALSE, lwd=2)
plot(ROC_rt, time=3, col='green', title=FALSE, lwd=2,add = T)
plot(ROC_rt, time=5, col='blue', title=FALSE, lwd=2,add = T)
legend('bottomright',
       c(paste0('AUC at 1 years: ',round(ROC_rt$AUC[1],2)),
         paste0('AUC at 3 years: ',round(ROC_rt$AUC[2],2)),
         paste0('AUC at 5 years: ',round(ROC_rt$AUC[3],2))),
       col=c('red','green','blue'),lwd=2,bty = 'n')
dev.off()

dat_km<-data.frame(time=GSE31210_module$time,
                   status=GSE31210_module$status,
                   riskscore=risk.GSE31210,
                   risk=ifelse(risk.GSE31210z>=0,'High','Low'))
write.table(dat_km,'geo_risk_module.txt',sep='\t',quote = F,row.names = F)

fit <- survfit(Surv(time, status) ~ risk, data = dat_km)

library(survminer)
pdf('geo_km.pdf',he=7,wi=7)
ggsurvplot(fit,risk.table=TRUE,
           conf.int=TRUE,
           palette = ggsci::pal_npg("nrc")(10),
           pval=TRUE,
           pval.method=TRUE,xlab='Year')
dev.off()
dat=data.frame(time=GSE31210_module$time,
               status=GSE31210_module$status,
               riskscore=risk.GSE31210z)
plot_cox_both<-function(dat,cutoff=0){
  colnames(dat)=c("time","status","riskscore")
  library(ggplot2)
  dat=dat[order(dat$riskscore),]
  dt1=data.frame(V1=1:length(dat$riskscore),V2=dat$riskscore,
                 RiskType=ifelse(dat$riskscore>cutoff,'High','Low')) 
  p1=ggplot(dt1, aes(x = V1, y = V2, colour = RiskType,fill=RiskType)) +
    geom_bar(stat = 'identity', position = 'dodge')+ggsci::scale_fill_npg()+theme_bw()+
    ylab('RiskScore')+
    theme(axis.text.y=element_text(family="Times",face="plain"),axis.text.x=element_blank()
          ,axis.title.x=element_blank(),legend.position='right'
          ,axis.ticks.x = element_blank()
          ,legend.background = element_rect(fill = NA, colour = NA)
          ,legend.title = element_text(family="Times",face="plain")
          ,legend.text = element_text(family="Times",face="plain"))
  dt2=data.frame(V1=1:length(dat$riskscore),V2=dat$time,
                 Status=ifelse(dat$status==1,'Dead','Alive'))  
  p2=ggplot(dt2, aes(x = V1, y = V2, colour = Status,shape =Status))+
    geom_point()+ggsci::scale_fill_npg()+theme_bw()+ylab('Time')+
    theme(axis.text.y=element_text(family="Times",face="plain"),
          axis.text.x=element_blank()
          ,axis.title.x=element_blank()
          ,legend.position='right'
          ,axis.ticks.x = element_blank()
          ,legend.background = element_rect(fill = NA, colour = NA)
          ,legend.title = element_text(family="Times",face="plain")
          ,legend.text = element_text(family="Times",face="plain"))
  p=ggarrange(p2,p1,ncol = 1,nrow = 2,align = "v")
  return(p)
}
plot_both<-plot_cox_both(dat,cutoff = 0)
plot_both
ggsave(filename = 'plot_both.pdf',plot = plot_both,he=7,wi=5)
