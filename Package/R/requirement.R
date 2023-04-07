#' requirement
#' @return Install required packages
#' @export
#'
#' @examples exvar::requirement()

requirement <- function(){
  # Detect operating system
  switch(Sys.info()[['sysname']],
         Windows= {sys <- "windows"},
         Linux  = {sys <- "linux"},
         Darwin = {sys <- "mac"}) # is Darwin correct ?
  # specie
  cat(paste0("These are the species currently supported by Exvar, choose the number corresponding to the target specie: \n",
             "[1] Homo sapiens \n",
             "[2] Mus musculus \n",
             "[3] Arabidopsis thaliana \n",
             "[4] Drosophila melanogaster \n",
             "[5] Danio rerio \n",
             "[6] Rattus norvegicus \n",
             "[7] Saccharomyces cerevisiae \n",
             "[8] Caenorhabditis elegans \n"))
  species <- readline("Type the number of the species that you would like to use as a reference: ")
  # install packages depending on operating system
  if (sys <-"linux"){
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
  switch(species,
         "1"={
           if (!require("BSgenome.Hsapiens.UCSC.hg19")) BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
           if (!require("BSgenome.Hsapiens.UCSC.hg19")) BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
           if (!require("org.Hs.eg.db")) BiocManager::install("org.Hs.eg.db")
           if (!require("TxDb.Hsapiens.UCSC.hg19.knownGene")) BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
           if (!require("SNPlocs.Hsapiens.dbSNP144.GRCh37")) BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37")
           if (!require("XtraSNPlocs.Hsapiens.dbSNP144.GRCh37")) BiocManager::install("XtraSNPlocs.Hsapiens.dbSNP144.GRCh37")
           if (!require("SNPlocs.Hsapiens.dbSNP151.GRCh38")) BiocManager::install("SNPlocs.Hsapiens.dbSNP151.GRCh38")
           if (!require("XtraSNPlocs.Hsapiens.dbSNP144.GRCh38")) BiocManager::install("XtraSNPlocs.Hsapiens.dbSNP144.GRCh38")
           
         },
         "2"={
           if (!require("BSgenome.Mmusculus.UCSC.mm10")) BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
           if (!require("org.Mm.eg.db")) BiocManager::install("org.Mm.eg.db")
           if (!require("TxDb.Mmusculus.UCSC.mm10.knownGene")) BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
           
         },
         "3"={
           if (!require("BSgenome.Athaliana.TAIR.TAIR9")) BiocManager::install("BSgenome.Athaliana.TAIR.TAIR9")
           if (!require("org.At.tair.db")) BiocManager::install("org.At.tair.db")
           if (!require("TxDb.Athaliana.BioMart.plantsmart28")) BiocManager::install("TxDb.Athaliana.BioMart.plantsmart28")
           
         },
         "4"={
           if (!require("BSgenome.Dmelanogaster.UCSC.dm6")) BiocManager::install("BSgenome.Dmelanogaster.UCSC.dm6")
           if (!require("org.Dm.eg.db")) BiocManager::install("org.Dm.eg.db")
           if (!require("TxDb.Dmelanogaster.UCSC.dm6.ensGene")) BiocManager::install("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
           
         },
         "5"={
           if (!require("BSgenome.Drerio.UCSC.danRer11")) BiocManager::install("BSgenome.Drerio.UCSC.danRer11")
           if (!require("org.Dr.eg.db")) BiocManager::install("org.Dr.eg.db")
           if (!require("TxDb.Drerio.UCSC.danRer11.refGene")) BiocManager::install("TxDb.Drerio.UCSC.danRer11.refGene")
           
         },
         "6"={ # rat
           if (!require("BSgenome.Rnorvegicus.UCSC.rn5")) BiocManager::install("BSgenome.Rnorvegicus.UCSC.rn5")
           if (!require("org.Rn.eg.db")) BiocManager::install("org.Rn.eg.db")
           if (!require("TxDb.Rnorvegicus.UCSC.rn5.refGene"))BiocManager::install("TxDb.Rnorvegicus.UCSC.rn5.refGene")
           
         },
         "7"={

           if (!require("BSgenome.Scerevisiae.UCSC.sacCer3")) BiocManager::install("BSgenome.Scerevisiae.UCSC.sacCer3")
           if (!require("org.Sc.sgd.db")) BiocManager::install("org.Sc.sgd.db")
           if (!require("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene")) BiocManager::install("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene")
           
         },
         "8"={
           if (!require("BSgenome.Celegans.UCSC.ce11")) BiocManager::install("BSgenome.Celegans.UCSC.ce11")
           if (!require("org.Ce.eg.db")) BiocManager::install("org.Ce.eg.db")
           if (!require("TxDb.Celegans.UCSC.ce11.refGene")) BiocManager::install("TxDb.Celegans.UCSC.ce11.refGene")

         }
  ) 
  } else {
    # put packages independent of OS
  }
  
}
