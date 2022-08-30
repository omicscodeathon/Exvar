.libPaths( c("/cbio/projects/003/imraan/r_lib", .libPaths()))
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rfastp")
BiocManager::install("ShortRead")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install("GenomeInfoDb")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
BiocManager::install("Rsubread")
BiocManager::install("Rsamtools")
BiocManager::install("GenomicAlignments")
BiocManager::install("SummarizedExperiment")
BiocManager::install("tools")
BiocManager::install("VariantTools")
BiocManager::install("VariantAnnotation")
BiocManager::install("saasCNV")
BiocManager::install("gmapR")
BiocManager::install("R.utils")
BiocManager::install("BiocParallel")

library("tcltk")
library("ShortRead")
library("Rfastp")
library("BSgenome.Hsapiens.UCSC.hg19")
library("GenomeInfoDb")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("Rsubread")
library("Rsamtools")
library("GenomicAlignments")
library("ggplot2")
library("SummarizedExperiment")
library("tools")
library("VariantTools")
library("VariantAnnotation")
library("saasCNV")
library("gmapR")
library("R.utils")
library("BiocParallel")

## This creates a fasta file containing the entire reference genome sans alt
## Resulting index is used for alignment
print("Creating reference genome index...")
refgen <- GmapGenome(BSgenome.Hsapiens.UCSC.hg19, create = TRUE, 
                     directory = "/cbio/projects/003/imraan")
print("Obtaining GSNAP parameters...")
snapParam <- GsnapParam(refgen, unique_only = TRUE,
                    molecule = "DNA", nthreads = 4)

## QC reads of multiple fastq files using rfastp 
## Files can be in .fastq or .gz format
## Then determines the number of QC iterations required to QC all selected files
## Each iteration creates a new folder to place output files correspondingly
inputFastq <- 
  "/cbio/projects/003/melissa/ICGNMD_WES_data/IC_SGS_00061/IC-SGS_00061_1.fastq.gz"
print("Receiving fastq file...")
fastqn <- length(inputFastq)
inputfl <- 1:fastqn

for (x in inputfl) {
  print("Quality checking fastq...")
  fastqPath <- file.path(inputFastq[c(x)])
  dir.create(paste0("/cbio/projects/003/imraan", 
                    "/", 
                    file_path_sans_ext(basename(inputFastq[c(x)]))))
  setwd(paste0("/cbio/projects/003/imraan", 
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

  sessionInfo()
## Alignment using the subjunc function
## Iterates over all previously selected files
  print("Aligning reads...")
  fastqPath <- file.path(inputFastq[c(x)])
  setwd(paste0("/cbio/projects/003/imraan", 
               "/", 
               file_path_sans_ext(basename(inputFastq[c(x)]))))
  print(getwd())
  print("Unzipping fastq.gz")
  gunzip(paste0(file_path_sans_ext(basename(inputFastq[c(x)])),
			                  '_quality_checked_R1.fastq.gz'))
  read <- paste0(file_path_sans_ext(basename(inputFastq[c(x)])), 
                  '_quality_checked_R1.fastq')
  print("GSNAP...")
  output <- gsnap(read, input_b = NULL, params = snapParam)
  print("Creating bam file...")
  bamfl <- as(output, "BamFile") 
  gc()

## Gene Counting iterating over all the bam files created previously
  print("Counting genes...")
  fastqPath <- file.path(inputFastq[c(x)])
  setwd(paste0("/cbio/projects/003/imraan", 
               "/", 
               file_path_sans_ext(basename(inputFastq[c(x)]))))
  geneExons <- exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene,by="gene")
  geneBam <- paste0(file_path_sans_ext(basename(inputFastq[c(x)])), 
                    "_quality_checked_R1.sam.bam")
  myBam <- BamFile(geneBam, 
                   yieldSize = 10000)
  GeneCounts <- summarizeOverlaps(geneExons,myBam,
                                  ignore.strand = TRUE)
  myGeneGR <- rowRanges(GeneCounts)
  GeneCounts_hg19 <- GeneCounts[all(seqnames(myGeneGR)),]
  GeneCountsMatrix <- assay(GeneCounts_hg19)+1
  Sum <- SummarizedExperiment(GeneCountsMatrix)
  myCounts <- data.frame(Counts =GeneCountsMatrix[,1])
  write.csv(mycounts, "Gene Counts.csv")

## Variant calling via VariantTools
  print("Calling variants...")
  parallelparam <- MulticoreParam(workers = 4)
  genome(repeats) <- genome(refgen)
  tallyParam <- TallyVariantsParam(refgen, mask = repeats, BPPARAM = parallelparam)
  variantcall <- callVariants(bamfl, tallyParam)
  print("Variant calling completed...")
  sampleNames(variantcall) <- 
    paste0(file_path_sans_ext(basename(inputFastq[c(x)])))
  mcols(variantcall) <- NULL
  vcf <- asVCF(variantcall)
  print("Writing vcf...")
  writeVcf(vcf, paste0(file_path_sans_ext(basename(inputFastq[c(x)])), ".vcf"), 
                       index = TRUE)
  
#CNV calling
  print("Calling copy number variants...")
  variant_df <- 
    vcf2txt(paste0(file_path_sans_ext(basename(inputFastq[c(x)])), ".vcf"))
  cnv_data <- cnv.data(variant_df)
  cnv_call <- cnv.call(cnv_data)
  cnvcf <- asVCF(cnv_call, 
                 paste0(file_path_sans_ext(basename(inputFastq[c(x)]))))
  writevcf(cnvcf, paste0(file_path_sans_ext(basename(inputFastq[c(x)])), 
                         ".vcf"),
           index = TRUE)
  
}



