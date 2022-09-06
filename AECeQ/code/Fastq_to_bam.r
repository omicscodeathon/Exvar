library("tcltk")
library("Rfastp")
library("BSgenome.Hsapiens.UCSC.hg19")
library("GenomeInfoDb")
library("GenomicAlignments")
library("tools")
library("gmapR")
library("R.utils")


## This creates a fasta file containing the entire reference genome
## Resulting index is used for alignment
## If index already exists in working directory, the script will skip to
## GsnapParam
print("Creating reference genome index...")
refgen <- GmapGenome(BSgenome.Hsapiens.UCSC.hg19, create = TRUE, 
                     directory = getwd())
print("Obtaining GSNAP parameters...")
snapParam <- GsnapParam(refgen, unique_only = TRUE,
                    molecule = "RNA", nthreads = 4)

## QC reads of multiple fastq files using rfastp 
## Files can be in .fastq or .gz format
## A new folder for outputs are created for each fastq file in the working dir
inputFastq <- tk_choose.files()
print("Receiving fastq file...")
fastqn <- length(inputFastq)
inputfl <- 1:fastqn

for (x in inputfl) {
  dir.create(paste0(getwd(), 
                    "/", 
                    file_path_sans_ext(basename(inputFastq[c(x)]))))
  wd <- getwd()
}

for (x in inputfl) {
  print("Quality checking fastq...")
  fastqPath <- file.path(inputFastq[c(x)])
  setwd(paste0(wd, 
               "/", 
               file_path_sans_ext(basename(inputFastq[c(x)]))))
  json_report <- rfastp(read1 = inputFastq[c(x)], 
                        outputFastq = 
                          paste0(file_path_sans_ext(basename(inputFastq[c(x)])), 
                                 '_quality_checked'),
                        thread = 4)
  QC <- qcSummary(json_report)
  write.csv(QC, "QC_summary.csv")
  gc()

## Alignment using the gsnap function
## Iterates over all previously selected files
## Output is an indexed bam file
  fastqPath <- file.path(inputFastq[c(x)])
  print("Unzipping fastq.gz")
  gunzip(paste0(file_path_sans_ext(basename(inputFastq[c(x)])),
			                  '_quality_checked_R1.fastq.gz'))
  read <- paste0(file_path_sans_ext(basename(inputFastq[c(x)])), 
                  '_quality_checked_R1.fastq')
  print("Aligning reads...")
  output <- gsnap(read, input_b = NULL, params = snapParam)
  print("Creating bam file...")
  bamfl <- as(output, "BamFile") 
  gc()

}



