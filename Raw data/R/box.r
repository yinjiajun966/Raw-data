setwd
library('ggplot2')
library('RColorBrewer')
library(ggpubr)
data<-read.csv("box_EMT.csv",header=TRUE)
a<-ggplot(data,aes(x=group,y=EMT,fill=group))+ #??fill=
  stat_boxplot(geom = "errorbar",width=0.15,aes(color="black"))+ 
  geom_boxplot(size=1.0,fill=c("#ffb6c1", "#87cefa"),outlier.fill="white",outlier.color="white")+ #size
  geom_jitter(aes(fill=group),width =0.2,shape = 21,size=1.5)
  scale_fill_manual(values = c("#ff4040", "#0000ff","#F0E442"))+
  scale_color_manual(values=c("black","black","black"))
  ggtitle(" ")
  theme_bw()+
  theme(legend.position="none", 
        axis.text.x=element_text(colour="black",family="Times",size=12),
        axis.text.y=element_text(family="Times",size=12,face="plain"),
        axis.title.y=element_text(family="Times",size = 12,face="plain"),        axis.title.x=element_text(family="Times",size = 12,face="plain"),
        plot.title = element_text(family="Times",size=12,face="bold",hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  stat_compare_means()+
  ylab("EMT value")+xlab("STAD")
pdf('EMT.pdf',width = 4,height = 8)
a
dev.off()