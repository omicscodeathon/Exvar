#' Convert a FASTQ file to an indexed BAM
#' 
#' This function takes in FASTQ files and performs quality control before
#' aligning to a reference genome. It assumes paired-end samples are of the same
#' file name with an underscore (_) and a number to signify different reads of
#' the same sample. Each sample's outputs will be stored in a separate
#' directory.
#'  
#' @param file A list of paths to FASTQ files. If no paths are entered, defaults to all fastq files in dir.
#' @param dir Output directory.
#' @param threads The number of cores to use in the process.
#' @param molecule A character string indicating either DNA or RNA samples.
#' @return A list of file paths to created BAM files
#' @export
processfastq <- function(file = list_files_with_exts(dir = dir, 
                                                     exts = "fastq"),
                            dir = getwd(), 
                            paired = FALSE, 
                            threads = 4L,
                            molecule = "RNA") {
  cat(paste0("These are the species currently supported by ExpVar: \n",
             "[1] Homo sapiens (hg19) \n",
             "[2] Homo sapiens (hg38) \n", 
             "[3] Mus musculus \n",
             "[4] Arabidopsis thaliana \n",
             "[5] Drosophila melanogaster \n",
             "[6] Danio rerio \n",
             "[7] Rattus norvegicus \n",
             "[8] Saccharomyces cerevisiae \n",
             "[9] Caenorhabditis elegans \n"))
  species <- readline("Type the number of the species that you would like to use as a reference: ")
  
  ##Sets the reference genome that corresponds to the species chosen by the user
  switch(species,
         "1"={
           ##Homo sapiens hg19
           organism <- BSgenome.Hsapiens.UCSC.hg19
         
           ##Selects hg19 as the reference genome
           ##If reference doesn't exist within package directory, create one
           if (dir.exists(paste0(find.package("ExpVar"), "/hg19"))) {
             refgen <- GmapGenome(organism,
                                  directory = find.package("ExpVar"))
           } else {
             print("Reference genome not found. Creating reference. This might take a while...")
             refgen <- GmapGenome(organism,
                                  directory = find.package("ExpVar"),
                                  create = TRUE)
           }
         },
         "2"={
           ##Homo sapiens hg38
           organism <- BSgenome.Hsapiens.UCSC.hg38
           
           ##Selects hg38 as the reference genome
           ##If reference doesn't exist within package directory, create one
           if (dir.exists(paste0(find.package("ExpVar"), "/hg38"))) {
             refgen <- GmapGenome(organism,
                                  directory = find.package("ExpVar"))
           } else {
             print("Reference genome not found. Creating reference. This might take a while...")
             refgen <- GmapGenome(organism,
                                  directory = find.package("ExpVar"),
                                  create = TRUE)
           }
         },
         "3"={
           ##Mus musculus mm10
           organism <- BSgenome.Mmusculus.UCSC.mm10
         
           ##Selects mm10 as the reference genome
           ##If reference doesn't exist within package directory, create one
           if (dir.exists(paste0(find.package("ExpVar"), "/mm10"))) {
             refgen <- GmapGenome(organism,
                                  directory = find.package("ExpVar"))
           } else {
             print("Reference genome not found. Creating reference. This might take a while...")
             refgen <- GmapGenome(organism,
                                  directory = find.package("ExpVar"),
                                  create = TRUE)
           }
         },
         "4"={
           ##Arabidopsis thaliana TAIR9
           organism <- BSgenome.Athaliana.TAIR.TAIR9
         
           ##Selects hg19 as the reference genome
           ##If reference doesn't exist within package directory, create one
           if (dir.exists(paste0(find.package("ExpVar"), "/TAIR9"))) {
             refgen <- GmapGenome(organism,
                                  directory = find.package("ExpVar"))
           } else {
             print("Reference genome not found. Creating reference. This might take a while...")
             refgen <- GmapGenome(organism,
                                  directory = find.package("ExpVar"),
                                  create = TRUE)
           }
         },
         "5"={
           ##Drosophilia melanogaster dm6
           organism <- BSgenome.Dmelanogaster.UCSC.dm6
         
           ##Selects dm6 as the reference genome
           ##If reference doesn't exist within package directory, create one
           if (dir.exists(paste0(find.package("ExpVar"), "/dm6"))) {
             refgen <- GmapGenome(organism,
                                  directory = find.package("ExpVar"))
           } else {
             print("Reference genome not found. Creating reference. This might take a while...")
             refgen <- GmapGenome(organism,
                                  directory = find.package("ExpVar"),
                                  create = TRUE)
           }
         },
         "6"={
           ##Danio rerio danRer11
           organism <- BSgenome.Drerio.UCSC.danRer11
         
           ##Selects danRer11 as the reference genome
           ##If reference doesn't exist within package directory, create one
           if (dir.exists(paste0(find.package("ExpVar"), "/danRer11"))) {
             refgen <- GmapGenome(organism,
                                  directory = find.package("ExpVar"))
           } else {
             print("Reference genome not found. Creating reference. This might take a while...")
             refgen <- GmapGenome(organism,
                                  directory = find.package("ExpVar"),
                                  create = TRUE)
           }
         },
         "7"={
           ##Rattus norvegicus rn5
           organism <- BSgenome.Rnorvegicus.UCSC.rn5
         
           ##Selects danRer11 as the reference genome
           ##If reference doesn't exist within package directory, create one
           if (dir.exists(paste0(find.package("ExpVar"), "/rn5"))) {
             refgen <- GmapGenome(organism,
                                  directory = find.package("ExpVar"))
           } else {
             print("Reference genome not found. Creating reference. This might take a while...")
             refgen <- GmapGenome(organism,
                                  directory = find.package("ExpVar"),
                                  create = TRUE)
           }
         },
         "8"={
           ##Saccharomyces cerevisiae sacCer3
           organism <- BSgenome.Scerevisiae.UCSC.sacCer3
         
           ##Selects sacCer3 as the reference genome
           ##If reference doesn't exist within package directory, create one
           if (dir.exists(paste0(find.package("ExpVar"), "/sacCer3"))) {
             refgen <- GmapGenome(organism,
                                  directory = find.package("ExpVar"))
           } else {
             print("Reference genome not found. Creating reference. This might take a while...")
             refgen <- GmapGenome(organism,
                                  directory = find.package("ExpVar"),
                                  create = TRUE)
           }
         },
         "9"={
           ##Caenorhabditis elagans
           organism <- BSgenome.Celegans.UCSC.ce11
         
           ##Selects ce11 as the reference genome
           ##If reference doesn't exist within package directory, create one
           if (dir.exists(paste0(find.package("ExpVar"), "/ce11"))) {
             refgen <- GmapGenome(organism,
                                  directory = find.package("ExpVar"))
           } else {
             print("Reference genome not found. Creating reference. This might take a while...")
             refgen <- GmapGenome(organism,
                                  directory = find.package("ExpVar"),
                                  create = TRUE)
           }
         }
)
      
  print("Obtaining GSNAP parameters...")
  snapParam <- GsnapParam(refgen, unique_only = TRUE,
                          molecule = molecule, nthreads = threads)
  bams <- c()
  
  if (isTRUE(paired)) {
    inputFastq <- file
    print("Receiving fastq file...")
    inputFastq <- sort(inputFastq)
    fastqn <- length(inputFastq)/2
    inputfl <- 1:fastqn
    foldernames <- gsub("\\_.*","", file_path_sans_ext(basename(inputFastq)))
    
    for (x in inputfl) {
      dir.create(paste0(dir, 
                        "/", 
                        foldernames[c(2*x)]))
    }
    
    for (x in inputfl) {
      print("Quality checking fastq...")
      fastqPath <- file.path(c(inputFastq[c(2*x-1)], inputFastq[c(2*x)]))
      read1 <- file_path_as_absolute(fastqPath[1])
      read2 <- file_path_as_absolute(fastqPath[2])
      setwd(paste0(dir, 
                   "/", 
                   foldernames[c(2*x)]))
      json_report <- rfastp(read1 = read1, read2 = read2,  
                            outputFastq = 
                              paste0(foldernames[2*x], 
                                     '_quality_checked'),
                            thread = threads, maxReadLength = 200L)
      QC <- qcSummary(json_report)
      write.csv(QC, "QC_summary.csv")
      gc()
      
      ## Alignment using the gsnap function
      ## Iterates over all previously selected files
      ## Output is an indexed bam file
      read1 <- paste0(foldernames[2*x-1], 
                      '_quality_checked_R1.fastq.gz')
      read2 <- paste0(foldernames[2*x], 
                      '_quality_checked_R2.fastq.gz')
      print("Aligning reads...")
      output <- gsnap(read1, read2, params = snapParam,
                      output = paste0(getwd(), foldernames[x]))
      print("Creating bam file...")
      bamfl <- as(output, "BamFile")
      append(bams, bamfl)
      gc()
    }
  } else {
    inputFastq <- file
    print("Receiving fastq file...")
    inputFastq <- sort(inputFastq)
    fastqn <- length(inputFastq)
    inputfl <- 1:fastqn
    foldernames <- gsub("\\_.*","", file_path_sans_ext(basename(inputFastq)))
    
    for (x in inputfl) {
      dir.create(paste0(dir, 
                        "/", 
                        foldernames[c(x)]))
    }
    
    for (x in inputfl) {
      print("Quality checking fastq...")
      fastqPath <- file.path(inputFastq[c(x)])
      read1 <- file_path_as_absolute(fastqPath)
      setwd(paste0(dir, 
                   "/", 
                   foldernames[c(x)]))
      json_report <- rfastp(read1 = read1, 
                            outputFastq = 
                              paste0(foldernames[x], 
                                     '_quality_checked'),
                            thread = threads)
      QC <- qcSummary(json_report)
      write.csv(QC, "QC_summary.csv")
      gc()
      
      ## Alignment using the gsnap function
      ## Iterates over all previously selected files
      ## Output is an indexed bam file
      print("Unzipping fastq.gz")
      gunzip(paste0(foldernames[x], 
                    '_quality_checked_R1.fastq.gz'))
      read <- paste0(foldernames[x], 
                     '_quality_checked_R1.fastq')
      print("Aligning reads...")
      output <- gsnap(read, input_b = NULL, params = snapParam,
                      output = paste0(getwd(), foldernames[x]))
      print("Creating bam file...")
      bamfl <- as(output, "BamFile") 
      append(bams, bamfl)
      gc()
    }
  }
  setwd(wd)
  return(bams)
}
           
