
#' Call copy number variants
#' 
#' This function calls copy number variants from sample BAM files compared to
#' control BAM files. It assumes that BAM files are stored in separate folders
#' as is created by fastqProcession(). This function requires that control BAM
#' files are provided. Once complete, it creates a CSV file containing copy
#' number information.
#' 
#' @param controldir The parent directory of the sample directories.
#' @param control The names of the folders in which control BAM files are. If NULL, all folders in controldir will be checked for BAM files.
#' @param experimentdir The parent directory of sample on which to investigate copy numbers.
#' @param experiment The names of the folders in which sample BAM files are. If NULL, all folders in experimentdir will be checked for BAM files.
#' @param bed A character string indicating BED file path or a TxDb object from which to extract a BED file.
#' @param outputdir The directory in which to place copy number call.
#' @return A data frame containing copy number calls.
#' @export
callcnv <- function(controldir,
                    experimentdir,
                    outputdir = getwd()) {
  if(Sys.info()[['sysname']] != "Linux"){
    message("This function is only available on Linux.")
    stop()
  }
  cat(paste0("These are the species currently supported by Exvar: \n",
             "[1] Homo sapiens (hg19) \n",
             "[2] Homo sapiens (hg38) \n"))
  species <- readline("Type the number of the species that you would like to use as a reference: ")
  
  wd <- getwd()
  ##Sets the reference genome that corresponds to the species chosen by the user
  switch(species,
         "1"={
           ##Homo sapiens hg19
           library(TxDb.Hsapiens.UCSC.hg19.knownGene)
           bed <- TxDb.Hsapiens.UCSC.hg19.knownGene
           
           bedname <- paste0(metadata(bed)[metadata(bed)$name == "Genome", 2], ".bed")
           if (file.exists(bedname)){
             countWindows <- panelcn.mops::getWindows(bedname)
           } else {
             outcome <- readline("BED file does not exist. Create one? [y/n]")
             if (outcome == "y") {
               rtracklayer::export.bed(bed, bedname)
               countWindows <- panelcn.mops::getWindows(bedname)
             } else {
               stop("BED file required to continue.")
             }
           }
         },
         "2"={
           ##Homo sapiens hg38
           library(TxDb.Hsapiens.UCSC.hg38.knownGene)
           bed <- TxDb.Hsapiens.UCSC.hg38.knownGene
           
           bedname <- paste0(metadata(bed)[metadata(bed)$name == "Genome", 2], ".bed")
           if (file.exists(bedname)){
             countWindows <- panelcn.mops::getWindows(bedname)
           } else {
             outcome <- readline("BED file does not exist. Create one? [y/n]")
             if (outcome == "y") {
               rtracklayer::export.bed(bed, bedname)
               countWindows <- panelcn.mops::getWindows(bedname)
             } else {
               stop("BED file required to continue.")
             }
           }
         }
  )
  
  ##Labels the folder that contains the sample folders under the path parameter
  control_samples <- list.dirs(path = controldir,
                               full.names = TRUE, recursive = FALSE)
  experiment_samples <- list.dirs(path = experimentdir, 
                                  full.names = TRUE, recursive = FALSE)
  
  ##Gets the BAM files from sample directories
  control_bamfl <- c()
  experiment_bamfl <- c()
  for (i in control_samples) {
    bam <- tools::list_files_with_exts(dir = i, exts = "bam")
    control_bamfl <- append(control_bamfl, bam)
  }
  
  for (i in experiment_samples) {
    bam <- tools::list_files_with_exts(dir = i, exts = "bam")
    experiment_bamfl <- append(experiment_bamfl, bam)
  }
  
  ##Counts copy number variation of regions found in the bed file
  ##The results are stored in a dataframe list with each index corresponding to
  ##each sample
  ##Control samples are excluded
  print("Getting gene counts...")
  a <- Rsamtools::BamFile(experiment_bamfl[1])
  countWindows <- countWindows[countWindows$chromosome == seqlevels(a),]
  print("Control")
  control <- panelcn.mops::countBamListInGRanges(countWindows = countWindows,
                                   bam.files = control_bamfl)
  #control <- getReadCountsFromBAM(control_bamfl, basename(control_bamfl), seqlevels(a), parallel = 10, min.mapq = 20, read.width = 200)
  print("Experiment")
  experiment <- panelcn.mops::countBamListInGRanges(countWindows = countWindows,
                                      bam.files = experiment_bamfl)
  #experiment <- getReadCountsFromBAM(experiment_bamfl, basename(experiment_bamfl), seqlevels(a), parallel = 10, min.mapq = 20, read.width = 200)
  index <- length(colnames(GenomicRanges::elementMetadata(experiment)))
  GenomicRanges::elementMetadata(experiment) <- cbind(GenomicRanges::elementMetadata(experiment),
                                       GenomicRanges::elementMetadata(control))
  print("Comparing copy numbers...")
  resultlist <- panelcn.mops::runPanelcnMops(experiment, 1:index,
                               countWindows = countWindows)
  sampleNames <- colnames(GenomicRanges::elementMetadata(experiment))
  resulttable <- panelcn.mops::createResultTable(resultlist = resultlist, XandCB = experiment,
                                   countWindows = countWindows,
                                   sampleNames = sampleNames)
  
  ##Creates a CSV file containing the copy number counts of all the samples
  ##CNVRanger has NE_id instead of sample_ID which is more generic
  print("Formatting results...")
  CNV_calls <- data.frame()
  for (data in resulttable) {
    data$Sample <- gsub("\\..*","", data$Sample)
    data$Sample <- gsub("_quality_checked_.*","", data$Sample)
    data$CN <- gsub("CN", "", data$CN)
    calls <- data.frame(chr = paste0("chr", data$Chr),
                        start = data$Start,
                        end = data$End,
                        sample_ID = data$Sample,
                        state = data$CN)
    CNV_calls <- merge(CNV_calls, calls, all.x = TRUE, all.y = TRUE)
    CNV_calls <- merge(CNV_calls, calls, all.x = TRUE, all.y = TRUE)
  }
  
  setwd(outputdir)
  print("Writing copy number variation to file...")
  write.csv(CNV_calls, "Copy Number Calls.csv", row.names = FALSE)
  CNV <- read.csv("Copy Number Calls.csv")
  setwd(wd)
  return(CNV)
}
