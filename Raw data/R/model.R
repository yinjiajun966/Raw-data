setwd
tcga_dat_m<-read.delim('tcga_dat_m.txt',sep='\t',row.names = 1,header = T,check.names = F)
tcga_dat_m[1:4,1:4]
library(survival)
tcga.cox=c()
for (i in colnames(tcga_dat_m)[3:ncol(tcga_dat_m)]){
  dat=data.frame(time=tcga_dat_m$OS.time,
                 status=tcga_dat_m$OS,
                 gene=tcga_dat_m[,i])
  fmla <- as.formula("Surv(time, status) ~gene")
  cox <- survival::coxph(fmla, data = dat)
  re = c(summary(cox)[[7]][5], 
         summary(cox)[[7]][2], 
         summary(cox)[[8]][3], 
         summary(cox)[[8]][4])
  tcga.cox=rbind(tcga.cox,re)
}
head(tcga.cox)
rownames(tcga.cox)=colnames(tcga_dat_m)[3:ncol(tcga_dat_m)]
colnames(tcga.cox) = c("p.value", "HR", "Low 95%CI", "High 95%CI")
tcga.cox=as.data.frame(tcga.cox)
write.csv(tcga.cox,"cox.csv")
sig.names=rownames(tcga.cox[which(tcga.cox$p.value < 0.05),])
length(sig.names)

library(glmnet)
set.seed(1111)

fit1=glmnet(as.matrix(tcga_dat_m[,sig.names])
            ,cbind(time=tcga_dat_m$OS.time,
                   status=tcga_dat_m$OS)
            ,family="cox"
            ,nlambda=100
            , alpha=1) 

cv.fit<-cv.glmnet(as.matrix(tcga_dat_m[,sig.names])
                  ,cbind(time=tcga_dat_m$OS.time,
                         status=tcga_dat_m$OS)
                  ,family="cox"
                  ,nlambda=100
                  , alpha=1)
sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
cv.fit$lambda.min
length(names(sig.coef))

pdf('lasso.pdf',width = 10,height = 5)
par(mfrow=c(1,2))
plot(fit1, xvar="lambda")
plot(cv.fit)
dev.off()


tcga_dat1 <- data.frame(time=tcga_dat_m$OS.time,
                        status=tcga_dat_m$OS,
                        tcga_dat_m[,names(sig.coef)])

fmla <- as.formula(paste0("Surv(time, status) ~"
                          ,paste0(names(sig.coef),collapse = '+')))

cox <- coxph(fmla, data =as.data.frame(tcga_dat1))
cox1 <- step(cox, trace = 0)
lan <- coef(cox1)
round(lan, 3)
genes <- names(cox1$coefficients)
length(genes)
paste0(round(lan, 3), '*', names(lan),collapse = '+')

risk.tcga <- as.numeric(lan%*%as.matrix(t(tcga_dat1[,genes])))

risk.tcgaz <- mosaic::zscore(risk.tcga)
tcga_dat1$time=tcga_dat1$time

write.table(genes,'gene.txt',sep='\t',quote = F,row.names = F,col.names = F)

library(pROC)
ROC_rt=timeROC::timeROC(T=tcga_dat1$time, 
                        delta=tcga_dat1$status,
                        marker=risk.tcgaz, cause=1,
                        weighting='marginal',
                        times=c(1,3,5), 
                        ROC=TRUE,iid = T)
ROC_rt$AUC
pdf('TCGA_ROC.pdf',he=7,wi=7)
plot(ROC_rt, time=1, col='red', title=FALSE, lwd=2)
plot(ROC_rt, time=3, col='green', title=FALSE, lwd=2,add = T)
plot(ROC_rt, time=5, col='blue', title=FALSE, lwd=2,add = T)
legend('bottomright',
       c(paste0('AUC at 1 years: ',round(ROC_rt$AUC[1],2)),
         paste0('AUC at 3 years: ',round(ROC_rt$AUC[2],2)),
         paste0('AUC at 5 years: ',round(ROC_rt$AUC[3],2))),
       col=c('red','green','blue'),lwd=2,bty = 'n')
dev.off()

dat_km<-data.frame(time=tcga_dat1$time,
                   status=tcga_dat1$status,
                   riskscore=risk.tcga,
                                      risk=ifelse(risk.tcgaz>=0,'High','Low'))
rownames(dat_km)=rownames(tcga_dat1)

write.table(dat_km,'tcga_risk_module.txt',sep='\t',quote = F,row.names = T)

fit <- survfit(Surv(time, status) ~ risk, data = dat_km)

library(survminer)
pdf('tcga_km.pdf',he=7,wi=7)
ggsurvplot(fit,risk.table=TRUE,
           conf.int=TRUE,
           palette = ggsci::pal_npg("nrc")(10),
           pval=TRUE,
           pval.method=TRUE,xlab='Year')
dev.off()
#
dat=data.frame(time=tcga_dat1$time,
               status=tcga_dat1$status,
               riskscore=risk.tcgaz)

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
a<-paste0(round(lan, 3), '*', names(lan),collapse = '+')
write.csv(a,"model.csv")





