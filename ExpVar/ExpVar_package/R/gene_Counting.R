#' Count genes
#' 
#' This function counts reads of gene regions between sample groups.
#' It assumes that sample BAM files are ordered in a directory structure such as
#' "group/sample/" as fastqProcession() would order it. It outputs a CSV file
#' showing gene counts. Works similarly to geneExpression(), but outputs count
#' data instead of differential expression data.
#' 
#' @param dir The parent directory of the sample groups. 
#' @param groups Folder names of the sample groups. The default is all folders in dir.
#' @param TxDb A TxDb object upon which regions of the genome are counted.
#' @param orgDb An orgDb object for annotating the CSV with gene symbols and Ensembl IDs.
#' @param outputdir Output directory of CSV file.
#' @param threads Number of cores to use.
#' @param paired Indicates whether the samples are from paired-end reads.
#' @return A data frame containing gene counts.
#' @export
gene_Counting <- function(dir = getwd(),
                           groups,
                           TxDb,
                           orgDb,
                           outputdir = getwd(),
                           threads = 4L,
                           paired = FALSE) {
  wd <- getwd()
  bpp = MulticoreParam(threads)
  geneExons <- exonsBy(TxDb, by = "gene")
  if (isTRUE(groups)) {
    groups <- paste0(dir, "/", groups)
    folders <- normalizePath(groups)
  } else {
    folders <- list.dirs(path = dir, full.names = TRUE, recursive = FALSE)
  }
  totaldir <- 1:length(folders)
  bamFilesToCount <- c()
  groupNames <- c()
  samples <- c()
  
  for (x in totaldir) {
    sampledir <- list.dirs(path = folders[c(x)], full.names = TRUE)
    sampledir <- sampledir[-1]
    totsample <- 1:length(sampledir)
    
    for (i in totsample) {
      print(basename(dirname(sampledir[c(i)])))
      sampletype <- basename(dirname(sampledir[c(i)]))
      groupNames <- append(groupNames, sampletype)
      print(tail(groupNames, n = 1L))
      bam <- list_files_with_exts(dir = sampledir[c(i)], exts = "bam")
      bamFilesToCount <- append(bamFilesToCount, bam)
      names <- basename(dirname(bam))
      samples <- append(samples, names)
    }
    
  }
  ##This creates a list of bam files to do gene counting with
  ##The folder where the bam files are found (ie sample folder) will inform the 
  ##name of the sample
  ##Samples are grouped according to the condition (named by the condition
  ##folders)
  
  ##A DESeq object is created from the gene counts
  ##If bam files aren't indexed, they are sorted and indexed first
  for (x in bamFilesToCount) {
    if (isTRUE(file.exists(paste0(x, ".bai")))) {
      setwd(dirname(x))
      sortBam(x, x)
      indexBam(x)
    }
  }
  setwd(wd)
  names(bamFilesToCount) <- samples
  myBams <- BamFileList(bamFilesToCount, yieldSize = 10000, asMates = paired)
  geneCounts <- summarizeOverlaps(geneExons, myBams, ignore.strand = TRUE,
                                  BPPARAM = bpp, singleEnd = isFALSE(paired))
  metaData <- data.frame(Group = groupNames, 
                         row.names = colnames(geneCounts))
  countMatrix <- assay(geneCounts)
  countDF <- data.frame(countMatrix)
  annotatedCount <- countDF
  
  if (isTRUE(orgDb)) {
    eToSym <- AnnotationDbi::select(orgDb,
                                    keys = rownames(countDF),
                                    keytype = "ENTREZID",
                                    columns= c("SYMBOL", "ENSEMBL"))
    annotatedCount <- merge(eToSym,countDF,
                            by.x=1,
                            by.y=0,
                            all.x=FALSE,
                            all.y=TRUE)
    annotatedCount <- annotatedCount[order(as.integer(annotatedCount$ENTREZID)),]
  }
  
  
  write.csv(annotatedCount, "Count Data.csv")
  setwd(wd)
  return(annotatedCount)
}


