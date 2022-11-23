#' Analyse differential gene expression
#' 
#' This function analyses differentially expressed genes between sample groups.
#' It assumes that sample BAM files are ordered in a directory structure such as
#' "group/sample/" as fastqProcession() would order it. There should be more
#' than one sample per group or else differential expression analysis won't work
#' It outputs a CSV file showing differential expression (ordered by p-value).
#' It works similarly to gene_Counting(), but then further analyses those counts
#' to obtain differential expression data.
#' 
#' 
#' @param dir The parent directory of the sample groups. 
#' @param groups Folder names of the sample groups. The default is all folders in dir.
#' @param TxDb A TxDb object upon which regions of the genome are counted.
#' @param orgDb An orgDb object for annotating the CSV with gene symbols and Ensembl IDs.
#' @param outputdir Output directory of CSV file.
#' @param threads Number of cores to use.
#' @param paired Indicates whether the samples are from paired-end reads.
#' @return A data frame list containing all of the differential expression comparison.
#' @export
geneExpression <- function(dir = getwd(),
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
  countGRanges <- rowRanges(geneCounts)
  ddse <- DESeqDataSetFromMatrix(countMatrix, colData = metaData, 
                                 design = ~Group, rowRanges = countGRanges)
  dds <- DESeq(ddse)
  
  normCounts <- counts(dds, normalized = TRUE)
  dispersion <- plotDispEsts(dds)
  ##Dispersion can be plotted with plotDispEsts(dds)
  
  ##Regardless of how many conditions there are, we can compare different
  ##conditions to each other as long as there are replicate samples
  groups <- unique(groupNames)
  groupmatrix <- combn(groups, 2)
  comparison <- length(groupmatrix)/2
  num_comparison <- 1:comparison
  compare <- c()
  lfc_compare <- c()
  DFcompare <- c()
  
  for (x in num_comparison) {
    myRes <- results(dds, contrast = c("Group", 
                                       groupmatrix[c(2*x)], 
                                       groupmatrix[c(2*x-1)]))
    myRes <- myRes[order(myRes$pvalue), ]
    compare <- append(compare, myRes)
    
    myResAsDF <- as.data.frame(myRes)
    myResAsDF$newPadj <- p.adjust(myResAsDF$pvalue)
    myResAsDF <- myResAsDF[!is.na(myResAsDF$padj), ]
    myResAsDF <- myResAsDF[order(myResAsDF$pvalue), ]
    annotatedRes <- myResAsDF
    
    if (isTRUE(orgDb)) {
      eToSym <- AnnotationDbi::select(orgDb,
                                      keys = rownames(myResAsDF),
                                      keytype = "ENTREZID",
                                      columns= c("SYMBOL", "ENSEMBL"))
      annotatedRes <- merge(eToSym,myResAsDF,
                            by.x=1,
                            by.y=0,
                            all.x=FALSE,
                            all.y=TRUE)
    }
    annotatedRes <- annotatedRes[order(annotatedRes$pvalue),]
    DFcompare <- append(DFcompare, annotatedRes)
    
    myRes_lfc <- lfcShrink(dds, coef = paste0("Group_",
                                              groupmatrix[c(2*x)],
                                              "_vs_",
                                              groupmatrix[c(2*x-1)]))
    lfc_compare <- append(lfc_compare, myRes_lfc)
    
    setwd(outputdir)
    write.csv(annotatedRes, paste0("Group_",
                                   groupmatrix[c(2*x)],
                                   "_vs_",
                                   groupmatrix[c(2*x-1)],
                                   ".csv"))
    
  }
  setwd(wd)
  return(DFcompare)
}

