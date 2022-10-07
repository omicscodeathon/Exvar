# dependencies

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require(c( "shiny","shinydashboard","DT","waiter"), quietly = TRUE))
  install.packages(c( 'shiny','shinydashboard','DT','waiter' ))

library(shiny)
library(shinydashboard)
library(waiter)
library(DT)


if (!require(c("CNVRanger","BSgenome.Hsapiens.UCSC.hg19","readr","AnnotationHub","Gviz","regioneR","BSgenome.Btaurus.UCSC.bosTau6.masked"), quietly = TRUE))
  BiocManager::install(c('AnnotationHub','readr','BSgenome.Hsapiens.UCSC.hg19','CNVRanger','Gviz','regioneR','BSgenome.Btaurus.UCSC.bosTau6.masked'))
  
  library(readr)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(CNVRanger)
  library(Gviz)
  library(AnnotationHub)
  library(regioneR)
  library(BSgenome.Btaurus.UCSC.bosTau6.masked)
