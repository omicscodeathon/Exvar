library(DESeq2)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BiocParallel)
library(Rsamtools)
library(GenomicFeatures)
library(GenomicAlignments)
library(org.Hs.eg.db)

##After creating bam files, each condition should have separate folders for 
##each of the samples
##This code checks all of these folders for the bam files contained within
##It also assumes each folder has a bam file
##Set working directory to the folder containing the different conditions folder
bpp = MulticoreParam(16)
geneExons <- exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by = "gene")
folders <- list.dirs(full.names = TRUE, recursive = FALSE)
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
names(bamFilesToCount) <- samples
myBams <- BamFileList(bamFilesToCount, yieldSize = 10000)
geneCounts <- summarizeOverlaps(geneExons, myBams, ignore.strand = TRUE,
                                BPPARAM = bpp)
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
  
  eToSym <- select(org.Hs.eg.db,
                   keys = rownames(myResAsDF),
                   keytype = "ENTREZID",
                   columns="SYMBOL")
  annotatedRes <- merge(eToSym,myResAsDF,
                        by.x=1,
                        by.y=0,
                        all.x=FALSE,
                        all.y=TRUE)
  annotatedRes <- annotatedRes[order(annotatedRes$pvalue),]
  DFcompare <- append(DFcompare, annotatedRes)
  
  myRes_lfc <- lfcShrink(dds, coef = paste0("Group_",
                                            groupmatrix[c(2*x)],
                                            "_vs_",
                                            groupmatrix[c(2*x-1)]))
  lfc_compare <- append(lfc_compare, myRes_lfc)
}
##Outputs contrast results as a dataframe
##Each result is indexed in DFcompare
##Genes are ordered according to p-value
##compare and lfc compare can be used to extract plots using plotMA