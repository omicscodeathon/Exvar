library("tcltk")
library("ShortRead")
library("Rfastp")
library("BSgenome.Hsapiens.UCSC.hg38")
library("GenomeInfoDb")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("Rsubread")
library("Rsamtools")
library("GenomicAlignments")
library("ggplot2")
library("SummarizedExperiment")
library("tools")
​
## This creates a fasta file containing the entire reference genome sans alt
## Resulting index is used for alignment
chr <- 1:22
chrxy <- c("X", "Y")
chr <- append(chr, chrxy)
mainChromosomes <- paste0("chr", chr)
mainChrSeq <- lapply(mainChromosomes,
                     function(x)BSgenome.Hsapiens.UCSC.hg38[[x]])
names(mainChrSeq) <- mainChromosomes
mainChrSeqSet <- DNAStringSet(mainChrSeq)
writeXStringSet(mainChrSeqSet,
                paste0("Hsapiens_hg38.fasta"))
gc()
myExons <- exons(TxDb.Hsapiens.UCSC.hg38.knownGene,
                 columns = c("tx_id", "exon_id"))
myExons <- myExons[seqnames(myExons) == mainChromosomes]
myExons <- myExons[lengths(myExons$exon_id) == 1]
dfExons <- as.data.frame(myExons)
print("Creating SAF")
dfExons$exon_id
SAF <- data.frame(GeneID = dfExons$exon_id,
                  Chr = dfExons$seqnames,
                  Start = dfExons$start,
                  End = dfExons$end,
                  Strand = dfExons$strand)
remove(myExons)
remove(dfExons)
gc()
dir.create(paste0(getwd(),"/Hsapiens_hg38"))
setwd(getwd(),"/Hsapiens_hg38")
buildindex(paste0("Hsapiens_hg38"),
           paste0("Hsapiens_hg38.fasta"), gappedIndex = TRUE, 
           indexSplit =TRUE, memory = 2000)
indexDir <- paste0(getwd(), "/Hsapiens_hg38")
gc()
​
## QC reads of multiple fastq files using rfastp 
## Files can be in .fastq or .gz format
## Then determines the number of QC iterations required to QC all selected files
## Each iteration creates a new folder to place output files correspondingly
inputFastq <- tk_choose.files(caption = "Select fastq files", 
                              filters = matrix(c(
                                "gzip", ".gz", 
                                "Fastq", ".fq",
                                "Fastq", ".fastq"), 3, 2, byrow = TRUE))
fastqn <- length(inputFastq)
inputfl <- 1:fastqn
​
for (x in inputfl) {
  fastqPath <- file.path(inputFastq[c(x)])
  dir.create(paste0(dirname(fastqPath), 
                    "/", 
                    file_path_sans_ext(basename(inputFastq[c(x)]))))
  setwd(paste0(dirname(fastqPath), 
               "/", 
               file_path_sans_ext(basename(inputFastq[c(x)]))))
  json_report <- rfastp(read1 = inputFastq[c(x)], 
                        outputFastq = 
                          paste0(file_path_sans_ext(basename(inputFastq[c(x)])), 
                                 '_quality_checked'),
                        thread = 4)
  QC <- qcSummary(json_report)
  write.table(QC, "QC_summary.txt", 
              sep = "\t", 
              row.names = TRUE, 
              col.names = TRUE)
  QCplot <-curvePlot(json_report)
  tiff(file = "QC_summary.tiff",
       width=6, height=4, units="in", res=100)
  plot(QCplot)
  dev.off()
  gc()
}
​
​
## Alignment using the subjunc function
## Iterates over all previously selected files
for (x in inputfl) {
  fastqPath <- file.path(inputFastq[c(x)])
  setwd(paste0(dirname(fastqPath), 
               "/", 
               file_path_sans_ext(basename(inputFastq[c(x)]))))
  print(getwd())
  read <- (paste0(file_path_sans_ext(basename(inputFastq[c(x)])), 
                  '_quality_checked_R1.fastq.gz'))
  outputBam <- paste0(file_path_sans_ext(basename(inputFastq[c(x)])), 
                      "_aligned.bam")
  subjunc(indexDir,
          read,
          output_file = 
            outputBam, 
          annot.ext = SAF, isGTF = FALSE,
          useAnnotation = TRUE)
  base <- file_path_sans_ext(basename(inputFastq[c(x)]))
  sortread <- paste0(base, "_aligned.bam")
  sortout <- paste0(base, "_sorted")
  sortBam(sortread, sortout)
  indexBam(paste0(sortout, ".bam"))
  gc()
}
​
## Gene Counting iterating over all the bam files created previously
for (x in inputfl) {
  fastqPath <- file.path(inputFastq[c(x)])
  setwd(paste0(dirname(fastqPath), 
               "/", 
               file_path_sans_ext(basename(inputFastq[c(x)]))))
  geneExons <- exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene,by="gene")
  geneBam <- paste0(file_path_sans_ext(basename(inputFastq[c(x)])), 
                    "_sorted.bam")
  myBam <- BamFile(geneBam, 
                   yieldSize = 10000)
  GeneCounts <- summarizeOverlaps(geneExons,myBam,
                                  ignore.strand = TRUE)
  myGeneGR <- rowRanges(GeneCounts)
  GeneCounts_hg38 <- GeneCounts[all(seqnames(myGeneGR) ==
                                      mainChromosomes),]
  GeneCountsMatrix <- assay(GeneCounts_hg38)+1
  Sum <- SummarizedExperiment(GeneCountsMatrix)
  Sum
  myCounts <- data.frame(Counts =GeneCountsMatrix[,1])
  tiff(file = "Gene Counts.tiff",
       width=6, height=4, units="in", res=100)
  ggplot(myCounts,aes(x=Counts)
  )+geom_density(fill="Red")+scale_x_log10()+theme_minimal() 
  dev.off()
}