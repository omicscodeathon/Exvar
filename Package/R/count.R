#' Count genes
#' 
#' This function counts reads of gene regions between sample groups.
#' It assumes that sample BAM files are ordered in a directory structure such as
#' "group/sample/" as processfastq() would order it. It outputs a CSV file
#' showing gene counts. Works similarly to expression(), but outputs count
#' data instead of differential expression data.
#' 
#' @param dir The parent directory of the sample groups. 
#' @param groups Folder names of the sample groups. The default is all folders in dir.
#' @param outputdir Output directory of CSV file.
#' @param threads Number of cores to use.
#' @param paired Indicates whether the samples are from paired-end reads.
#' @export
#' @return A data frame containing gene counts.

count <- function(dir = getwd(),
                  groups = NULL,
                  outputdir = getwd(),
                  threads = 4L,
                  paired = FALSE) {
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
  
  wd <- getwd()
  bpp <- BiocParallel::MulticoreParam(threads)
  ##Sets the reference genome that corresponds to the species chosen by the user
  switch(species,
         "1"={
           ##Homo sapiens hg19
           library(BSgenome.Hsapiens.UCSC.hg19)
           library(TxDb.Hsapiens.UCSC.hg19.knownGene)
           library(org.Hs.eg.db)
           organism <- BSgenome.Hsapiens.UCSC.hg19
           
           geneExons <- GenomicFeatures::exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by = "gene")
           if (is.null(groups)) {
             folders <- list.dirs(path = dir, full.names = TRUE, recursive = FALSE)
           } else {
             groups <- paste0(dir, "/", groups)
             folders <- normalizePath(groups)
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
               bam <- tools::list_files_with_exts(dir = sampledir[c(i)], exts = "bam")
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
           names(bamFilesToCount) <- samples
           myBams <- Rsamtools::BamFileList(bamFilesToCount, yieldSize = 10000, asMates = paired)
           geneCounts <- GenomicAlignments::summarizeOverlaps(geneExons, myBams, ignore.strand = TRUE,
                                           BPPARAM = bpp, singleEnd = isFALSE(paired))
           metaData <- data.frame(Group = groupNames, 
                                  row.names = colnames(geneCounts))
           countMatrix <- SummarizedExperiment::assay(geneCounts)
           countDF <- data.frame(countMatrix)
           annotatedCount <- countDF
           
           if (isTRUE(orgDb)) {
             eToSym <- AnnotationDbi::select(org.Hs.eg.db,
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
         },
         "2"={
           ##Homo sapiens hg38
           library(BSgenome.Hsapiens.UCSC.hg38)
           library(TxDb.Hsapiens.UCSC.hg38.knownGene)
           library(org.Hs.eg.db)
           organism <- BSgenome.Hsapiens.UCSC.hg38
           
           geneExons <- GenomicFeatures::exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = "gene")
           if (is.null(groups)) {
             folders <- list.dirs(path = dir, full.names = TRUE, recursive = FALSE)
           } else {
             groups <- paste0(dir, "/", groups)
             folders <- normalizePath(groups)
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
               bam <- tools::list_files_with_exts(dir = sampledir[c(i)], exts = "bam")
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
           names(bamFilesToCount) <- samples
           myBams <- Rsamtools::BamFileList(bamFilesToCount, yieldSize = 10000, asMates = paired)
           geneCounts <- GenomicAlignments::summarizeOverlaps(geneExons, myBams, ignore.strand = TRUE,
                                           BPPARAM = bpp, singleEnd = isFALSE(paired))
           metaData <- data.frame(Group = groupNames, 
                                  row.names = colnames(geneCounts))
           countMatrix <- SummarizedExperiment::assay(geneCounts)
           countDF <- data.frame(countMatrix)
           annotatedCount <- countDF
           
           if (isTRUE(orgDb)) {
             eToSym <- AnnotationDbi::select(org.Hs.eg.db,
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
         },
         "3"={
           ##Mus musculus mm10
           library(BSgenome.Mmusculus.UCSC.mm10)
           library(TxDb.Mmusculus.UCSC.mm10.knownGene)
           library(org.Mm.eg.db)
           organism <- BSgenome.Mmusculus.UCSC.mm10
           
           geneExons <- GenomicFeatures::exonsBy(TxDb.Mmusculus.UCSC.mm10.knownGene, by = "gene")
           if (is.null(groups)) {
             folders <- list.dirs(path = dir, full.names = TRUE, recursive = FALSE)
           } else {
             groups <- paste0(dir, "/", groups)
             folders <- normalizePath(groups)
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
               bam <- tools::list_files_with_exts(dir = sampledir[c(i)], exts = "bam")
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
           names(bamFilesToCount) <- samples
           myBams <- Rsamtools::BamFileList(bamFilesToCount, yieldSize = 10000, asMates = paired)
           geneCounts <- GenomicAlignments::summarizeOverlaps(geneExons, myBams, ignore.strand = TRUE,
                                           BPPARAM = bpp, singleEnd = isFALSE(paired))
           metaData <- data.frame(Group = groupNames, 
                                  row.names = colnames(geneCounts))
           countMatrix <- SummarizedExperiment::assay(geneCounts)
           countDF <- data.frame(countMatrix)
           annotatedCount <- countDF
           
           if (isTRUE(orgDb)) {
             eToSym <- AnnotationDbi::select(org.Mm.eg.db,
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
         },
         "4"={
           ##Arabidopsis thaliana TAIR9
           library(BSgenome.Athaliana.TAIR.TAIR9)
           library(TxDb.Athaliana.BioMart.plantsmart28)
           library(org.At.tair.db)
           organism <- BSgenome.Athaliana.TAIR.TAIR9
           
           geneExons <- GenomicFeatures::exonsBy(xDb.Athaliana.BioMart.plantsmart28, by = "gene")
           if (is.null(groups)) {
             folders <- list.dirs(path = dir, full.names = TRUE, recursive = FALSE)
           } else {
             groups <- paste0(dir, "/", groups)
             folders <- normalizePath(groups)
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
               bam <- tools::list_files_with_exts(dir = sampledir[c(i)], exts = "bam")
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
           names(bamFilesToCount) <- samples
           myBams <- Rsamtools::BamFileList(bamFilesToCount, yieldSize = 10000, asMates = paired)
           geneCounts <- GenomicAlignments::summarizeOverlaps(geneExons, myBams, ignore.strand = TRUE,
                                           BPPARAM = bpp, singleEnd = isFALSE(paired))
           metaData <- data.frame(Group = groupNames, 
                                  row.names = colnames(geneCounts))
           countMatrix <- SummarizedExperiment::assay(geneCounts)
           countDF <- data.frame(countMatrix)
           annotatedCount <- countDF
           
           if (isTRUE(orgDb)) {
             eToSym <- AnnotationDbi::select(org.At.tair.db,
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
         },
         "5"={
           ##Drosophilia melanogaster dm6
           library(BSgenome.Dmelanogaster.UCSC.dm6)
           library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
           library(org.Dm.eg.db)
           organism <- BSgenome.Dmelanogaster.UCSC.dm6
           
           geneExons <- GenomicFeatures::exonsBy(TxDb.Dmelanogaster.UCSC.dm6.ensGene, by = "gene")
           if (is.null(groups)) {
             folders <- list.dirs(path = dir, full.names = TRUE, recursive = FALSE)
           } else {
             groups <- paste0(dir, "/", groups)
             folders <- normalizePath(groups)
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
               bam <- tools::list_files_with_exts(dir = sampledir[c(i)], exts = "bam")
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
           names(bamFilesToCount) <- samples
           myBams <- Rsamtools::BamFileList(bamFilesToCount, yieldSize = 10000, asMates = paired)
           geneCounts <- GenomicAlignments::summarizeOverlaps(geneExons, myBams, ignore.strand = TRUE,
                                           BPPARAM = bpp, singleEnd = isFALSE(paired))
           metaData <- data.frame(Group = groupNames, 
                                  row.names = colnames(geneCounts))
           countMatrix <- SummarizedExperiment::assay(geneCounts)
           countDF <- data.frame(countMatrix)
           annotatedCount <- countDF
           
           if (isTRUE(orgDb)) {
             eToSym <- AnnotationDbi::select(org.Dm.eg.db,
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
         },
         "6"={
           ##Danio rerio danRer11
           library(BSgenome.Drerio.UCSC.danRer11)
           library(TxDb.Drerio.UCSC.danRer11.refGene)
           library(org.Dr.eg.db)
           organism <- BSgenome.Drerio.UCSC.danRer11
           
           geneExons <- GenomicFeatures::exonsBy(TxDb.Drerio.UCSC,danRer11.refGene, by = "gene")
           if (is.null(groups)) {
             folders <- list.dirs(path = dir, full.names = TRUE, recursive = FALSE)
           } else {
             groups <- paste0(dir, "/", groups)
             folders <- normalizePath(groups)
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
               bam <- tools::list_files_with_exts(dir = sampledir[c(i)], exts = "bam")
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
           names(bamFilesToCount) <- samples
           myBams <- Rsamtools::BamFileList(bamFilesToCount, yieldSize = 10000, asMates = paired)
           geneCounts <- GenomicAlignments::summarizeOverlaps(geneExons, myBams, ignore.strand = TRUE,
                                           BPPARAM = bpp, singleEnd = isFALSE(paired))
           metaData <- data.frame(Group = groupNames, 
                                  row.names = colnames(geneCounts))
           countMatrix <- SummarizedExperiment::assay(geneCounts)
           countDF <- data.frame(countMatrix)
           annotatedCount <- countDF
           
           if (isTRUE(orgDb)) {
             eToSym <- AnnotationDbi::select(org.Dr.eg.db,
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
         },
         "7"={
           ##Rattus norvegicus rn5
           library(BSgenome.Rnorvegicus.UCSC.rn5)
           library(TxDb.Dnorvegicus.UCSC.rn5.refGene)
           library(org.Rn.eg.db)
           organism <- BSgenome.Rnorvegicus.UCSC.rn5
           
           geneExons <- GenomicFeatures::exonsBy(TxDb.Rnorvegicus.UCSC.rn5.refGene, by = "gene")
           if (is.null(groups)) {
             folders <- list.dirs(path = dir, full.names = TRUE, recursive = FALSE)
           } else {
             groups <- paste0(dir, "/", groups)
             folders <- normalizePath(groups)
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
               bam <- tools::list_files_with_exts(dir = sampledir[c(i)], exts = "bam")
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
           names(bamFilesToCount) <- samples
           myBams <- Rsamtools::BamFileList(bamFilesToCount, yieldSize = 10000, asMates = paired)
           geneCounts <- GenomicAlignments::summarizeOverlaps(geneExons, myBams, ignore.strand = TRUE,
                                           BPPARAM = bpp, singleEnd = isFALSE(paired))
           metaData <- data.frame(Group = groupNames, 
                                  row.names = colnames(geneCounts))
           countMatrix <- SummarizedExperiment::assay(geneCounts)
           countDF <- data.frame(countMatrix)
           annotatedCount <- countDF
           
           if (isTRUE(orgDb)) {
             eToSym <- AnnotationDbi::select(org.Rn.eg.db,
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
         },
         "8"={
           ##Saccharomyces cerevisiae sacCer3
           library(BSgenome.Scerevisiae.UCSC.sacCer3)
           library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
           library(org.Sc.sgd.db)
           organism <- BSgenome.Scerevisiae.UCSC.sacCer3
           
           geneExons <- GenomicFeatures::exonsBy(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene, by = "gene")
           if (is.null(groups)) {
             folders <- list.dirs(path = dir, full.names = TRUE, recursive = FALSE)
           } else {
             groups <- paste0(dir, "/", groups)
             folders <- normalizePath(groups)
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
               bam <- tools::list_files_with_exts(dir = sampledir[c(i)], exts = "bam")
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
           names(bamFilesToCount) <- samples
           myBams <- Rsamtools::BamFileList(bamFilesToCount, yieldSize = 10000, asMates = paired)
           geneCounts <- GenomicAlignments::summarizeOverlaps(geneExons, myBams, ignore.strand = TRUE,
                                           BPPARAM = bpp, singleEnd = isFALSE(paired))
           metaData <- data.frame(Group = groupNames, 
                                  row.names = colnames(geneCounts))
           countMatrix <- SummarizedExperiment::assay(geneCounts)
           countDF <- data.frame(countMatrix)
           annotatedCount <- countDF
           
           if (isTRUE(orgDb)) {
             eToSym <- AnnotationDbi::select(org.Sc.sgd.db,
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
         },
         "9"={
           ##Caenorhabditis elagans
           library(BSgenome.Celegans.UCSC.ce11)
           library(TxDb.Delegans.UCSC.ce11.refGene)
           library(org.Ce.eg.db)
           organism <- BSgenome.Celegans.UCSC.ce11
           
           geneExons <- GenomicFeatures::exonsBy(TxDb.Celegans.UCSC.ce11.refGene, by = "gene")
           if (is.null(groups)) {
             folders <- list.dirs(path = dir, full.names = TRUE, recursive = FALSE)
           } else {
             groups <- paste0(dir, "/", groups)
             folders <- normalizePath(groups)
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
               bam <- tools::list_files_with_exts(dir = sampledir[c(i)], exts = "bam")
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
           names(bamFilesToCount) <- samples
           myBams <- Rsamtools::BamFileList(bamFilesToCount, yieldSize = 10000, asMates = paired)
           geneCounts <- GenomicAlignments::summarizeOverlaps(geneExons, myBams, ignore.strand = TRUE,
                                           BPPARAM = bpp, singleEnd = isFALSE(paired))
           metaData <- data.frame(Group = groupNames, 
                                  row.names = colnames(geneCounts))
           countMatrix <- SummarizedExperiment::assay(geneCounts)
           countDF <- data.frame(countMatrix)
           annotatedCount <- countDF
           
           if (isTRUE(orgDb)) {
             eToSym <- AnnotationDbi::select(org.Ce.eg.db,
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
         }
  )
  write.csv(annotatedCount, "Count Data.csv")
  setwd(wd)
  return(annotatedCount)
}
