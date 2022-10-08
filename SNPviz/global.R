# dependencies

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

############  app 
#if (!require(c( "shiny","shinydashboard","DT","waiter"), quietly = TRUE))
#  install.packages(c( 'shiny','shinydashboard','DT','waiter' ))

library(shiny)
library(shinydashboard)
library(waiter)
library(DT)


#############   expression part
#if (!require(c( "dplyr","ggplot2","gplots","ggrepel", "echarts4r"), quietly = TRUE))
#  install.packages(c( 'dplyr','ggplot2','gplots','ggrepel', 'echarts4r'))

#if (!require(c("limma","DESeq2","AnnotationDbi","org.Mm.eg.db","org.hs.eg.db",
#               "ReportingTools","GO.db","GOstats","pathview","gage","gageData","select"), quietly = TRUE))
#  BiocManager::install(c('limma','DESeq2','AnnotationDbi','org.Mm.eg.db','org.hs.eg.db',
#                         'ReportingTools','GO.db','GOstats','pathview','gage','gageData','select'))
  
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
  library(ComplexHeatmap)
  #ontology
  library(ClusterProfiler)
  library(AnnotationDbi)
  library(org.Hs.eg.db)

############### SNP part 

#if (!require("VariantAnnotation", quietly = TRUE))
 # BiocManager::install("VariantAnnotation")
  library(echarts4r)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(GenomicFeatures)
  library(stringr)