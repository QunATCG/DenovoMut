{
  rm(list = ls())
  setwd("/Users/liqun/Desktop/")
  library(ggplot2)
  go = read.table("./Goplot.txt", header = T, sep = "\t")
  #go = go[go$ONTOLOGY == "BP",]
  go = go[go$Pvalue<=0.01,]
  go = go[order(go$Pvalue, decreasing=F),]
  go = go[1:5,]
  go = go[!duplicated(go$Description),]
  #go = go[c(-5,-8,-11),]
  head(go)
  
  ggplot(data=go, aes(x=reorder(Description,-Pvalue),y=-log10(Pvalue))) + 
    geom_bar(stat="identity", width=0.8, fill='#AEC28E') + 
    coord_flip() + 
    theme_bw() + theme(panel.grid=element_blank()) 
}