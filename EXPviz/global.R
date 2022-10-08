# dependencies

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

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
#BiocManager::install("clusterProfiler", force = TRUE)

library(readr)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
#ontology
#library(ClusterProfiler,lib.loc = "/usr/lib/R/library")
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)


