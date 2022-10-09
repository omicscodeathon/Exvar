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
#' @param genome A BSgenome object, GmapGenome object, or a character string indicating the genome name eg. "hg19".
#' @param genomedir A directory containing the reference genome. Otherwise, it is the parent directory of the reference genome is genome is a character string or BSgenome object.
#' @param paired Indicates whether the samples are from paired-end or single-end reads.
#' @param threads The number of cores to use in the process.
#' @param molecule A character string indicating either DNA or RNA samples.
#' @return A list of file paths to created BAM files
#' @export
fastqProcession <- function(file = list_files_with_exts(dir = dir, 
                                                        exts = "fastq"),
                            dir = getwd(), 
                            genome,
                            genomedir,
                            paired = FALSE, 
                            threads = 4L,
                            molecule = "RNA") {
  wd <- getwd()
  if (class(genome) != "GmapGenome") {
    if(isFALSE(exists("genome"))) {
      if (exists("genomedir")) {
        if (dir.exists(genomedir)) {
          refgen <- GmapGenome(BSgenome.Hsapiens.UCSC.hg19, 
                               directory = dirname(genomedir),
                               name = basename(genomedir))
        }
      }
    } else {
      if (isFALSE(exists("genomedir"))) {
        genomedir <- getwd()
      }
      if (class(genome) == "character") {
        if (dir.exists(genome)) {
          refgen <- GmapGenome(BSgenome.Hsapiens.UCSC.hg19, genomedir, 
                               name = genome)
        } 
      } else {
        if (dir.exists(paste0(genomdir, "/", metadata(genome)$genome))) {
          refgen <- GmapGenome(genome, genomedir, create = FALSE) 
        } else {
          outcome <- readline("Reference genome doesn't exist in directory. Create one? [y/n]")
          if (outcome == "y") {
            refgen <- GmapGenome(genome, genomedir, create = TRUE)
          } else {
            stop("Cannot proceed without reference genome.")
          }
        }
      }
    }
  } else {
    refgen <- genome
  }

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
      output <- gsnap(read1, read2, params = snapParam)
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
      output <- gsnap(read, input_b = NULL, params = snapParam)
      print("Creating bam file...")
      bamfl <- as(output, "BamFile") 
      appened(bams, bamfl)
      gc()
    }
  }
  setwd(wd)
  return(bams)
}
