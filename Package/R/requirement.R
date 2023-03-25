#' requirement
#'
#' @param specie specie name.
#'
#' @return Install required packages
#' @export
#'
#' @examples expvar::requirement("Homo Sapiens")

requirement <- function(specie){
  #for all species
  if (!require("shiny")) install.packages("shiny")
  if (!require("shinydashboard")) install.packages("shinydashboard")
  if (!require("DT")) install.packages("DT")
  if (!require("shinyWidgets")) install.packages("shinyWidgets")
  if (!require("shinythemes")) install.packages("shinythemes")
  if (!require("dplyr")) install.packages("dplyr")
  if (!require("BiocParallel")) BiocManager::install("BiocParallel")
  if (!require("GenomeInfoDb")) BiocManager::install("GenomeInfoDb")
  if (!require("GenomicAlignments")) BiocManager::install("GenomicAlignments")
  if (!require("gmapR")) BiocManager::install("gmapR")
  if (!require("panelcn.mops")) BiocManager::install("panelcn.mops")
  if (!require("R.utils")) install.packages("R.utils")
  if (!require("Rfastp")) BiocManager::install("Rfastp")
  if (!require("Rsamtools"))  BiocManager::install("Rsamtools")
  if (!require("VariantTools"))  BiocManager::install("VariantTools")
  if (!require("tools")) install.packages("tools") # ?
  if (!require("data.table")) install.packages("data.table")
  if (!require("readr")) install.packages("readr")
  if (!require("DESeq2")) BiocManager::install("DESeq2")
  if (!require("ggplot2")) install.packages("ggplot2")
  if (!require("ComplexHeatmap")) BiocManager::install("ComplexHeatmap")
  if (!require("clusterProfiler")) BiocManager::install("clusterProfiler")
  if (!require("AnnotationDbi")) BiocManager::install("AnnotationDbi")
  if (!require("enrichplot")) BiocManager::install("enrichplot")
  if (!require("echarts4r")) install.packages("echarts4r")
  if (!require("GenomicFeatures"))  BiocManager::install("GenomicFeatures")
  if (!require("stringr")) install.packages("stringr")
  if (!require("VariantAnnotation"))  BiocManager::install("VariantAnnotation")
  if (!require("CNVRanger")) BiocManager::install("CNVRanger")
  if (!require("Gviz")) BiocManager::instal("Gviz")
  if (!require("AnnotationHub")) BiocManager::instal("AnnotationHub")
  if (!require("regioneR")) BiocManager::instal("regioneR")
  if (!require("ggnewscale")) install.packages("ggnewscale")
  if (!require("tibble")) install.packages("tibble")
  #species specific
  if (species == "Homo Sapiens"){
    if (!require("BSgenome.Hsapiens.UCSC.hg19")) BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
    if (!require("BSgenome.Hsapiens.UCSC.hg19")) BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
    if (!require("org.Hs.eg.db")) BiocManager::install("org.Hs.eg.db")
    if (!require("TxDb.Hsapiens.UCSC.hg19.knownGene")) BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
    if (!require("SNPlocs.Hsapiens.dbSNP144.GRCh37")) BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37")
    if (!require("XtraSNPlocs.Hsapiens.dbSNP144.GRCh37")) BiocManager::install("XtraSNPlocs.Hsapiens.dbSNP144.GRCh37")
    if (!require("SNPlocs.Hsapiens.dbSNP151.GRCh38")) BiocManager::install("SNPlocs.Hsapiens.dbSNP151.GRCh38")
    if (!require("XtraSNPlocs.Hsapiens.dbSNP144.GRCh38")) BiocManager::install("XtraSNPlocs.Hsapiens.dbSNP144.GRCh38")
    
  } else if (species == "Mus Musculus"){
    if (!require("BSgenome.Mmusculus.UCSC.mm10")) BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
    if (!require("org.Mm.eg.db")) BiocManager::install("org.Mm.eg.db")
    if (!require("TxDb.Mmusculus.UCSC.mm10.knownGene")) BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")

    
  } else if (species == "Arabidopsis Thaliana"){
    if (!require("BSgenome.Athaliana.TAIR.TAIR9")) BiocManager::install("BSgenome.Athaliana.TAIR.TAIR9")
    if (!require("org.At.tair.db")) BiocManager::install("org.At.tair.db")
    if (!require("TxDb.Athaliana.BioMart.plantsmart28")) BiocManager::install("TxDb.Athaliana.BioMart.plantsmart28")

    
  }  else if (species == "Drosophila Melanogaster"){
    if (!require("BSgenome.Dmelanogaster.UCSC.dm6")) BiocManager::install("BSgenome.Dmelanogaster.UCSC.dm6")
    if (!require("org.Dm.eg.db")) BiocManager::install("org.Dm.eg.db")
    if (!require("TxDb.Dmelanogaster.UCSC.dm6.ensGene")) BiocManager::install("TxDb.Dmelanogaster.UCSC.dm6.ensGene")

    
  } else if (species == "Danio rerio"){
    if (!require("BSgenome.Drerio.UCSC.danRer11")) BiocManager::install("BSgenome.Drerio.UCSC.danRer11")
    if (!require("org.Dr.eg.db")) BiocManager::install("org.Dr.eg.db")
    if (!require("TxDb.Drerio.UCSC.danRer11.refGene")) BiocManager::install("TxDb.Drerio.UCSC.danRer11.refGene")

    
  } else if (species == "Rattus_norvegicus"){
    if (!require("BSgenome.Rnorvegicus.UCSC.rn5")) BiocManager::install("BSgenome.Rnorvegicus.UCSC.rn5")
    if (!require("org.Rn.eg.db")) BiocManager::install("org.Rn.eg.db")
    if (!require("TxDb.Rnorvegicus.UCSC.rn5.refGene"))BiocManager::install("TxDb.Rnorvegicus.UCSC.rn5.refGene")

    
  } else if (species == "Saccharomyces Cerevisiae"){
    if (!require("BSgenome.Scerevisiae.UCSC.sacCer3")) BiocManager::install("BSgenome.Scerevisiae.UCSC.sacCer3")
    if (!require("org.Sc.sgd.db")) BiocManager::install("org.Sc.sgd.db")
    if (!require("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene")) BiocManager::install("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene")

    
  } else {
    print("Non-Valid specie!")
  }
  
}
