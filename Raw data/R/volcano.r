setwd
library(tidyverse)
library(ggplot2)
library(ggrepel)
data <- read.csv("volcano.csv",row.names = 1)
logFC = 1
adj.P.Val = 0.01
k1 <- (data$adj.P.Val < adj.P.Val) & (data$logFC < -logFC)
k2 <- (data$adj.P.Val < adj.P.Val) & (data$logFC > logFC)
data <- mutate(data, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))
data$symbol <- rownames(data)
data$label <- NA  

p <- ggplot(data = data, 
            aes(x = logFC, 
                y = -log10(adj.P.Val))) +  
  geom_point(alpha = 0.4, size = 3.5, 
             aes(color = change)) +     
  ylab("-log10(adj.P.Val)") +              
  scale_color_manual(values = c("blue4", "grey", "red3")) + 
  geom_vline(xintercept = c(-logFC, logFC), lty = 4, col = "black", lwd = 0.8) + 
  geom_hline(yintercept = -log10(adj.P.Val), lty = 4, col = "black", lwd = 0.8) +   
  theme_bw() 
p

p0 <- p + geom_label_repel(data = data, aes(label = label),
                           size = 4,                          
                           box.padding = unit(0.5, "lines"),  
                           point.padding = unit(0.8, "lines"),
                           segment.color = "black",           
                           show.legend = FALSE,               
                           max.overlaps = 10000)              
p0
pdf('volcano.pdf',width = 8,height = 8)
p0
dev.off()