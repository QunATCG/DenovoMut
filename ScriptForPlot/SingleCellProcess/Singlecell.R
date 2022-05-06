## clear env
{
  rm(list = ls())
  setwd("/Users/liqun/Desktop/Qun/DenovoControl/")
  outputDir = getwd()
}
##
{
  cran.packages <- c("msigdbr", "dplyr", "purrr", "stringr","magrittr",
                     "RobustRankAggreg", "tibble", "reshape2", "ggsci",
                     "tidyr", "aplot", "ggfun", "ggplotify", "ggridges",
                     "gghalves", "Seurat", "SeuratObject", "methods",
                     "devtools", "BiocManager","data.table","doParallel",
                     "doRNG")
  if (!requireNamespace(cran.packages, quietly = TRUE)) {
    install.packages(cran.packages, ask = F, update = F)
  }
  
  # install packages from Bioconductor
  bioconductor.packages <- c("GSEABase", "AUCell", "SummarizedExperiment",
                             "singscore", "GSVA", "ComplexHeatmap", "ggtree",
                             "Nebulosa")
  if (!requireNamespace(bioconductor.packages, quietly = TRUE)) {
    BiocManager::install(bioconductor.packages, ask = F, update = F)
  }
}


## Libraries
{
  suppressMessages(library(Seurat))
  suppressMessages(library(ggplot2))
  suppressMessages(library(gridExtra))
  suppressMessages(library(dplyr))
  suppressMessages(library(plotly))
  suppressMessages(library(genefilter))
  suppressMessages(library(gplots))
  #BiocManager::install("tidyverse")
  #library(tidyverse)
}

# FGC.rds can download from https://zenodo.org/record/6523440/files/FGC.rds?download=1
FGC <- readRDS(file = "./FGC.rds")

gonadalSomaticMakerGenes <- c("WT1", "COL3A1", "COL1A1", "VCAN", "FOXL2",
                              "ALDH1A1")
immuneMarkerGenes <- c("CD68", "CD4")
endothelialMarkerGenes <- c("VWF", "CD34")
FGCMakerGenes <- c("KIT", "DAZL")
earlyFGCMarkerGenes <- c("GSTM1", "GABRA3", "POU5F1")
preMeioticLateFGCMarkerGenes <- c("STRA8")
meioticFGCMarkerGenes <- c("SPO11")
dictyateOocyteMarkerGenes <- c("ZP3")

geneSymbol = c(gonadalSomaticMakerGenes,immuneMarkerGenes, FGCMakerGenes, earlyFGCMarkerGenes, preMeioticLateFGCMarkerGenes, meioticFGCMarkerGenes,dictyateOocyteMarkerGenes)

# example plot
FeaturePlot(object = FGC, 
            features = dictyateOocyteMarkerGenes, 
            cols = c("grey", "blue"), 
            reduction = "tsne")
