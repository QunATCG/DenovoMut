setwd('C:/Users/liqun/Desktop/Qun/3_estimateRisk/Code/')
data = read.table('toPlot.txt',header = T, sep = " ")
library(ggplot2)

ggplot_settings=function(legend=F){
  pp=theme(plot.title = element_text(hjust = 0.5,size=rel(1),face="bold"),
           panel.border = element_blank(), 
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(), 
           panel.background = element_blank(),
           
           strip.text = element_text(size=rel(1),face="bold"),
           strip.background = element_blank(),
           
           axis.text.x =element_text(size=rel(0.8),angle=0),
           axis.text.y =element_text(size=rel(0.8),angle=0),
           axis.title=element_text(size=rel(0.8),face="bold"),
           axis.line = element_line(colour = "black", size = rel(0.7), linetype = "solid")
  )
  
  if(legend==F){
    pp=theme(legend.position = "none")+pp
  }
  return(pp)
}

pdf("riskgene.pdf",width = 4.33333,height = 3.177083)
ggplot(data, mapping = aes(x = nGenes, y = likelihood)) + geom_point() + ggplot_settings()+
  geom_vline(aes(xintercept=223), colour="#BB0000", linetype="dashed") + 
  annotate('text', color= "blue",x=400, y=0.0001, label='MLE=223 Genes') +
  xlab("nGenes") + ylab("Denisty") + ggtitle("Maximum likelihood estimate of risk gene")
dev.off()
