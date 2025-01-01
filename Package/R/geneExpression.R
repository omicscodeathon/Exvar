
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
#' @param outputdir Output directory of CSV file.
#' @param threads Number of cores to use.
#' @param paired Indicates whether the samples are from paired-end reads.
#' @return A data frame list containing all of the differential expression comparison.
#' @export
geneExpression <- function(dir = getwd(),
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
  
  switch(species,
         "1"={
           ## Homo sapiens hg19
           library(BSgenome.Hsapiens.UCSC.hg19)
           library(TxDb.Hsapiens.UCSC.hg19.knownGene)
           library(org.Hs.eg.db)
           organism <- BSgenome.Hsapiens.UCSC.hg19
           
           wd <- getwd()
           bpp = BiocParallel::MulticoreParam(threads)
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
           countGRanges <- SummarizedExperiment::rowRanges(geneCounts)
           ddse <- DESeq2::DESeqDataSetFromMatrix(countMatrix, colData = metaData, 
                                          design = ~Group, rowRanges = countGRanges)
           dds <- DESeq2::DESeq(ddse)
           
           normCounts <- DESeq2::counts(dds, normalized = TRUE)
           dispersion <- DESeq2::plotDispEsts(dds)
           ##Dispersion can be plotted with DESeq2::plotDispEsts(dds)
           
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
             myRes <- DESeq2::results(dds, contrast = c("Group", 
                                                groupmatrix[c(2*x)], 
                                                groupmatrix[c(2*x-1)]))
             myRes <- myRes[order(myRes$pvalue), ]
             compare <- append(compare, myRes)
             
             myResAsDF <- as.data.frame(myRes)
             myResAsDF$newPadj <- stats::p.adjust(myResAsDF$pvalue)
             myResAsDF <- myResAsDF[!is.na(myResAsDF$padj), ]
             myResAsDF <- myResAsDF[order(myResAsDF$pvalue), ]
             annotatedRes <- myResAsDF
             
             
             eToSym <- AnnotationDbi::select(org.Hs.eg.db,
                                             keys = rownames(myResAsDF),
                                             keytype = "ENTREZID",
                                             columns= c("SYMBOL", "ENSEMBL"))
             annotatedRes <- merge(eToSym,myResAsDF,
                                   by.x=1,
                                   by.y=0,
                                   all.x=FALSE,
                                   all.y=TRUE)
             
             annotatedRes <- annotatedRes[order(annotatedRes$SYMBOL, annotatedRes$padj),]
             annotatedRes <- annotatedRes[!duplicated(annotatedRes$SYMBOL),]
             annotatedRes <- annotatedRes[order(annotatedRes$padj),]
             DFcompare <- append(DFcompare, annotatedRes)
             
             myRes_lfc <- DESeq2::lfcShrink(dds, coef = paste0("Group_",
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
         },
         "2"={
           ## Homo sapiens hg38
           library(BSgenome.Hsapiens.UCSC.hg38)
           library(TxDb.Hsapiens.UCSC.hg38.knownGene)
           library(org.Hs.eg.db)
           organism <- BSgenome.Hsapiens.UCSC.hg38
           
           wd <- getwd()
           bpp = BiocParallel::MulticoreParam(threads)
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
           countGRanges <- SummarizedExperiment::rowRanges(geneCounts)
           ddse <- DESeq2::DESeqDataSetFromMatrix(countMatrix, colData = metaData, 
                                          design = ~Group, rowRanges = countGRanges)
           dds <- DESeq2::DESeq(ddse)
           
           normCounts <- DESeq2::counts(dds, normalized = TRUE)
           dispersion <- DESeq2::plotDispEsts(dds)
           ##Dispersion can be plotted with DESeq2::plotDispEsts(dds)
           
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
             myRes <- DESeq2::results(dds, contrast = c("Group", 
                                                groupmatrix[c(2*x)], 
                                                groupmatrix[c(2*x-1)]))
             myRes <- myRes[order(myRes$pvalue), ]
             compare <- append(compare, myRes)
             
             myResAsDF <- as.data.frame(myRes)
             myResAsDF$newPadj <- stats::p.adjust(myResAsDF$pvalue)
             myResAsDF <- myResAsDF[!is.na(myResAsDF$padj), ]
             myResAsDF <- myResAsDF[order(myResAsDF$pvalue), ]
             annotatedRes <- myResAsDF
             
             
             eToSym <- AnnotationDbi::select(org.Hs.eg.db,
                                             keys = rownames(myResAsDF),
                                             keytype = "ENTREZID",
                                             columns= c("SYMBOL", "ENSEMBL"))
             annotatedRes <- merge(eToSym,myResAsDF,
                                   by.x=1,
                                   by.y=0,
                                   all.x=FALSE,
                                   all.y=TRUE)
             
             annotatedRes <- annotatedRes[order(annotatedRes$SYMBOL, annotatedRes$padj),]
             annotatedRes <- annotatedRes[!duplicated(annotatedRes$SYMBOL),]
             annotatedRes <- annotatedRes[order(annotatedRes$padj),]
             DFcompare <- append(DFcompare, annotatedRes)
             
             myRes_lfc <- DESeq2::lfcShrink(dds, coef = paste0("Group_",
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
         },
         "3"={
           ##Mus musculus mm10
           library(BSgenome.Mmusculus.UCSC.mm10)
           library(TxDb.Mmusculus.UCSC.mm10.knownGene)
           library(org.Mm.eg.db)
           organism <- BSgenome.Mmusculus.UCSC.mm10
           
           wd <- getwd()
           bpp = BiocParallel::MulticoreParam(threads)
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
           countGRanges <- SummarizedExperiment::rowRanges(geneCounts)
           ddse <- DESeq2::DESeqDataSetFromMatrix(countMatrix, colData = metaData, 
                                          design = ~Group, rowRanges = countGRanges)
           dds <- DESeq2::DESeq(ddse)
           
           normCounts <- DESeq2::counts(dds, normalized = TRUE)
           dispersion <- DESeq2::plotDispEsts(dds)
           ##Dispersion can be plotted with DESeq2::plotDispEsts(dds)
           
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
             myRes <- DESeq2::results(dds, contrast = c("Group", 
                                                groupmatrix[c(2*x)], 
                                                groupmatrix[c(2*x-1)]))
             myRes <- myRes[order(myRes$pvalue), ]
             compare <- append(compare, myRes)
             
             myResAsDF <- as.data.frame(myRes)
             myResAsDF$newPadj <- stats::p.adjust(myResAsDF$pvalue)
             myResAsDF <- myResAsDF[!is.na(myResAsDF$padj), ]
             myResAsDF <- myResAsDF[order(myResAsDF$pvalue), ]
             annotatedRes <- myResAsDF
             
             
             eToSym <- AnnotationDbi::select(org.Mm.eg.db,
                                             keys = rownames(myResAsDF),
                                             keytype = "ENTREZID",
                                             columns= c("SYMBOL", "ENSEMBL"))
             annotatedRes <- merge(eToSym,myResAsDF,
                                   by.x=1,
                                   by.y=0,
                                   all.x=FALSE,
                                   all.y=TRUE)
             
             annotatedRes <- annotatedRes[order(annotatedRes$SYMBOL, annotatedRes$padj),]
             annotatedRes <- annotatedRes[!duplicated(annotatedRes$SYMBOL),]
             annotatedRes <- annotatedRes[order(annotatedRes$padj),]
             DFcompare <- append(DFcompare, annotatedRes)
             
             myRes_lfc <- DESeq2::lfcShrink(dds, coef = paste0("Group_",
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
         },
         "4"={
           ##Arabidopsis thaliana TAIR9
           library(BSgenome.Athaliana.TAIR.TAIR9)
           library(TxDb.Athaliana.BioMart.plantsmart28)
           library(org.At.tair.db)
           organism <- BSgenome.Athaliana.TAIR.TAIR9
           
           wd <- getwd()
           bpp = BiocParallel::MulticoreParam(threads)
           geneExons <- GenomicFeatures::exonsBy(TxDb.Athaliana.BioMart.plantsmart28, by = "gene")
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
           countGRanges <- SummarizedExperiment::rowRanges(geneCounts)
           ddse <- DESeq2::DESeqDataSetFromMatrix(countMatrix, colData = metaData, 
                                          design = ~Group, rowRanges = countGRanges)
           dds <- DESeq2::DESeq(ddse)
           
           normCounts <- DESeq2::counts(dds, normalized = TRUE)
           dispersion <- DESeq2::plotDispEsts(dds)
           ##Dispersion can be plotted with DESeq2::plotDispEsts(dds)
           
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
             myRes <- DESeq2::results(dds, contrast = c("Group", 
                                                groupmatrix[c(2*x)], 
                                                groupmatrix[c(2*x-1)]))
             myRes <- myRes[order(myRes$pvalue), ]
             compare <- append(compare, myRes)
             
             myResAsDF <- as.data.frame(myRes)
             myResAsDF$newPadj <- stats::p.adjust(myResAsDF$pvalue)
             myResAsDF <- myResAsDF[!is.na(myResAsDF$padj), ]
             myResAsDF <- myResAsDF[order(myResAsDF$pvalue), ]
             annotatedRes <- myResAsDF
             
             eToSym <- AnnotationDbi::select(org.At.tair.db,
                                             keys = rownames(myResAsDF),
                                             keytype = "ENTREZID",
                                             columns= c("SYMBOL", "ENSEMBL"))
             annotatedRes <- merge(eToSym,myResAsDF,
                                   by.x=1,
                                   by.y=0,
                                   all.x=FALSE,
                                   all.y=TRUE)
             
             annotatedRes <- annotatedRes[order(annotatedRes$SYMBOL, annotatedRes$padj),]
             annotatedRes <- annotatedRes[!duplicated(annotatedRes$SYMBOL),]
             annotatedRes <- annotatedRes[order(annotatedRes$padj),]
             DFcompare <- append(DFcompare, annotatedRes)
             
             myRes_lfc <- DESeq2::lfcShrink(dds, coef = paste0("Group_",
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
         },
         "5"={
           ##Drosophilia melanogaster dm6
           library(BSgenome.Dmelanogaster.UCSC.dm6)
           library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
           library(org.Dm.eg.db)
           organism <- BSgenome.Dmelanogaster.UCSC.dm6
           
           wd <- getwd()
           bpp = BiocParallel::MulticoreParam(threads)
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
           countGRanges <- SummarizedExperiment::rowRanges(geneCounts)
           ddse <- DESeq2::DESeqDataSetFromMatrix(countMatrix, colData = metaData, 
                                          design = ~Group, rowRanges = countGRanges)
           dds <- DESeq2::DESeq(ddse)
           
           normCounts <- DESeq2::counts(dds, normalized = TRUE)
           dispersion <- DESeq2::plotDispEsts(dds)
           ##Dispersion can be plotted with DESeq2::plotDispEsts(dds)
           
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
             myRes <- DESeq2::results(dds, contrast = c("Group", 
                                                groupmatrix[c(2*x)], 
                                                groupmatrix[c(2*x-1)]))
             myRes <- myRes[order(myRes$pvalue), ]
             compare <- append(compare, myRes)
             
             myResAsDF <- as.data.frame(myRes)
             myResAsDF$newPadj <- stats::p.adjust(myResAsDF$pvalue)
             myResAsDF <- myResAsDF[!is.na(myResAsDF$padj), ]
             myResAsDF <- myResAsDF[order(myResAsDF$pvalue), ]
             annotatedRes <- myResAsDF
             
             eToSym <- AnnotationDbi::select(org.Dm.eg.db,
                                             keys = rownames(myResAsDF),
                                             keytype = "ENTREZID",
                                             columns= c("SYMBOL", "ENSEMBL"))
             annotatedRes <- merge(eToSym,myResAsDF,
                                   by.x=1,
                                   by.y=0,
                                   all.x=FALSE,
                                   all.y=TRUE)
             
             annotatedRes <- annotatedRes[order(annotatedRes$SYMBOL, annotatedRes$padj),]
             annotatedRes <- annotatedRes[!duplicated(annotatedRes$SYMBOL),]
             annotatedRes <- annotatedRes[order(annotatedRes$padj),]
             DFcompare <- append(DFcompare, annotatedRes)
             
             myRes_lfc <- DESeq2::lfcShrink(dds, coef = paste0("Group_",
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
         },
         "6"={
           ##Danio rerio danRer11
           library(BSgenome.Drerio.UCSC.danRer11)
           library(TxDb.Drerio.UCSC.danRer11.refGene)
           library(org.Dr.eg.db)
           organism <- BSgenome.Drerio.UCSC.danRer11
           
           wd <- getwd()
           bpp = BiocParallel::MulticoreParam(threads)
           geneExons <- GenomicFeatures::exonsBy(TxDb.Drerio.UCSC.danRer11.refGene, by = "gene")
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
           countGRanges <- SummarizedExperiment::rowRanges(geneCounts)
           ddse <- DESeq2::DESeqDataSetFromMatrix(countMatrix, colData = metaData, 
                                          design = ~Group, rowRanges = countGRanges)
           dds <- DESeq2::DESeq(ddse)
           
           normCounts <- DESeq2::counts(dds, normalized = TRUE)
           dispersion <- DESeq2::plotDispEsts(dds)
           ##Dispersion can be plotted with DESeq2::plotDispEsts(dds)
           
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
             myRes <- DESeq2::results(dds, contrast = c("Group", 
                                                groupmatrix[c(2*x)], 
                                                groupmatrix[c(2*x-1)]))
             myRes <- myRes[order(myRes$pvalue), ]
             compare <- append(compare, myRes)
             
             myResAsDF <- as.data.frame(myRes)
             myResAsDF$newPadj <- stats::p.adjust(myResAsDF$pvalue)
             myResAsDF <- myResAsDF[!is.na(myResAsDF$padj), ]
             myResAsDF <- myResAsDF[order(myResAsDF$pvalue), ]
             annotatedRes <- myResAsDF
             
             eToSym <- AnnotationDbi::select(org.Dr.eg.db,
                                             keys = rownames(myResAsDF),
                                             keytype = "ENTREZID",
                                             columns= c("SYMBOL", "ENSEMBL"))
             annotatedRes <- merge(eToSym,myResAsDF,
                                   by.x=1,
                                   by.y=0,
                                   all.x=FALSE,
                                   all.y=TRUE)
             
             annotatedRes <- annotatedRes[order(annotatedRes$SYMBOL, annotatedRes$padj),]
             annotatedRes <- annotatedRes[!duplicated(annotatedRes$SYMBOL),]
             annotatedRes <- annotatedRes[order(annotatedRes$padj),]
             DFcompare <- append(DFcompare, annotatedRes)
             
             myRes_lfc <- DESeq2::lfcShrink(dds, coef = paste0("Group_",
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
         },
         "7"={
           ##Rattus norvegicus rn5
           library(BSgenome.Rnorvegicus.UCSC.rn5)
           library(TxDb.Dnorvegicus.UCSC.rn5.refGene)
           library(org.Rn.eg.db)
           organism <- BSgenome.Rnorvegicus.UCSC.rn5
           
           wd <- getwd()
           bpp = BiocParallel::MulticoreParam(threads)
           geneExons <- GenomicFeatures::exonsBy(TxDb.Dnorvegicus.UCSC.rn5.refGene, by = "gene")
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
           countGRanges <- SummarizedExperiment::rowRanges(geneCounts)
           ddse <- DESeq2::DESeqDataSetFromMatrix(countMatrix, colData = metaData, 
                                          design = ~Group, rowRanges = countGRanges)
           dds <- DESeq2::DESeq(ddse)
           
           normCounts <- DESeq2::counts(dds, normalized = TRUE)
           dispersion <- DESeq2::plotDispEsts(dds)
           ##Dispersion can be plotted with DESeq2::plotDispEsts(dds)
           
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
             myRes <- DESeq2::results(dds, contrast = c("Group", 
                                                groupmatrix[c(2*x)], 
                                                groupmatrix[c(2*x-1)]))
             myRes <- myRes[order(myRes$pvalue), ]
             compare <- append(compare, myRes)
             
             myResAsDF <- as.data.frame(myRes)
             myResAsDF$newPadj <- stats::p.adjust(myResAsDF$pvalue)
             myResAsDF <- myResAsDF[!is.na(myResAsDF$padj), ]
             myResAsDF <- myResAsDF[order(myResAsDF$pvalue), ]
             annotatedRes <- myResAsDF
             
             eToSym <- AnnotationDbi::select(org.Rn.eg.db,
                                             keys = rownames(myResAsDF),
                                             keytype = "ENTREZID",
                                             columns= c("SYMBOL", "ENSEMBL"))
             annotatedRes <- merge(eToSym,myResAsDF,
                                   by.x=1,
                                   by.y=0,
                                   all.x=FALSE,
                                   all.y=TRUE)
             
             annotatedRes <- annotatedRes[order(annotatedRes$SYMBOL, annotatedRes$padj),]
             annotatedRes <- annotatedRes[!duplicated(annotatedRes$SYMBOL),]
             annotatedRes <- annotatedRes[order(annotatedRes$padj),]
             DFcompare <- append(DFcompare, annotatedRes)
             
             myRes_lfc <- DESeq2::lfcShrink(dds, coef = paste0("Group_",
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
         },
         "8"={
           ##Saccharomyces cerevisiae sacCer3
           library(BSgenome.Scerevisiae.UCSC.sacCer3)
           library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
           library(org.Sc.sgd.db)
           organism <- BSgenome.Scerevisiae.UCSC.sacCer3
           
           wd <- getwd()
           bpp = BiocParallel::MulticoreParam(threads)
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
           countGRanges <- SummarizedExperiment::rowRanges(geneCounts)
           ddse <- DESeq2::DESeqDataSetFromMatrix(countMatrix, colData = metaData, 
                                          design = ~Group, rowRanges = countGRanges)
           dds <- DESeq2::DESeq(ddse)
           
           normCounts <- DESeq2::counts(dds, normalized = TRUE)
           dispersion <- DESeq2::plotDispEsts(dds)
           ##Dispersion can be plotted with DESeq2::plotDispEsts(dds)
           
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
             myRes <- DESeq2::results(dds, contrast = c("Group", 
                                                groupmatrix[c(2*x)], 
                                                groupmatrix[c(2*x-1)]))
             myRes <- myRes[order(myRes$pvalue), ]
             compare <- append(compare, myRes)
             
             myResAsDF <- as.data.frame(myRes)
             myResAsDF$newPadj <- stats::p.adjust(myResAsDF$pvalue)
             myResAsDF <- myResAsDF[!is.na(myResAsDF$padj), ]
             myResAsDF <- myResAsDF[order(myResAsDF$pvalue), ]
             annotatedRes <- myResAsDF
             
             eToSym <- AnnotationDbi::select(org.Sc.sgd.db,
                                             keys = rownames(myResAsDF),
                                             keytype = "ENTREZID",
                                             columns= c("SYMBOL", "ENSEMBL"))
             annotatedRes <- merge(eToSym,myResAsDF,
                                   by.x=1,
                                   by.y=0,
                                   all.x=FALSE,
                                   all.y=TRUE)
             
             annotatedRes <- annotatedRes[order(annotatedRes$SYMBOL, annotatedRes$padj),]
             annotatedRes <- annotatedRes[!duplicated(annotatedRes$SYMBOL),]
             annotatedRes <- annotatedRes[order(annotatedRes$padj),]
             DFcompare <- append(DFcompare, annotatedRes)
             
             myRes_lfc <- DESeq2::lfcShrink(dds, coef = paste0("Group_",
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
         },
         "9"={
           ##Caenorhabditis elagans
           library(BSgenome.Celegans.UCSC.ce11)
           library(TxDb.Delegans.UCSC.ce11.refGene)
           library(org.Ce.eg.db)
           organism <- BSgenome.Celegans.UCSC.ce11
           
           wd <- getwd()
           bpp = BiocParallel::MulticoreParam(threads)
           geneExons <- GenomicFeatures::exonsBy(TxDb.Delegans.UCSC.ce11.refGene, by = "gene")
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
           countGRanges <- SummarizedExperiment::rowRanges(geneCounts)
           ddse <- DESeq2::DESeqDataSetFromMatrix(countMatrix, colData = metaData, 
                                          design = ~Group, rowRanges = countGRanges)
           dds <- DESeq2::DESeq(ddse)
           
           normCounts <- DESeq2::counts(dds, normalized = TRUE)
           dispersion <- DESeq2::plotDispEsts(dds)
           ##Dispersion can be plotted with DESeq2::plotDispEsts(dds)
           
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
             myRes <- DESeq2::results(dds, contrast = c("Group", 
                                                groupmatrix[c(2*x)], 
                                                groupmatrix[c(2*x-1)]))
             myRes <- myRes[order(myRes$pvalue), ]
             compare <- append(compare, myRes)
             
             myResAsDF <- as.data.frame(myRes)
             myResAsDF$newPadj <- stats::p.adjust(myResAsDF$pvalue)
             myResAsDF <- myResAsDF[!is.na(myResAsDF$padj), ]
             myResAsDF <- myResAsDF[order(myResAsDF$pvalue), ]
             annotatedRes <- myResAsDF
             
             eToSym <- AnnotationDbi::select(org.Ce.eg.db,
                                             keys = rownames(myResAsDF),
                                             keytype = "ENTREZID",
                                             columns= c("SYMBOL", "ENSEMBL"))
             annotatedRes <- merge(eToSym,myResAsDF,
                                   by.x=1,
                                   by.y=0,
                                   all.x=FALSE,
                                   all.y=TRUE)
             
             annotatedRes <- annotatedRes[order(annotatedRes$SYMBOL, annotatedRes$padj),]
             annotatedRes <- annotatedRes[!duplicated(annotatedRes$SYMBOL),]
             annotatedRes <- annotatedRes[order(annotatedRes$padj),]
             DFcompare <- append(DFcompare, annotatedRes)
             
             myRes_lfc <- DESeq2::lfcShrink(dds, coef = paste0("Group_",
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
         }
  )
  setwd(wd)
  return(DFcompare)
}
