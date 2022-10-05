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

##A csv file containing the count data is created
names(bamFilesToCount) <- samples
myBams <- BamFileList(bamFilesToCount, yieldSize = 10000)
geneCounts <- summarizeOverlaps(geneExons, myBams, ignore.strand = TRUE,
                                BPPARAM = bpp)
metaData <- data.frame(Group = groupNames, 
                       row.names = colnames(geneCounts))
countMatrix <- assay(geneCounts)
countDF <- data.frame(countMatrix)
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

write.csv(annotatedCount, "Count Data.csv")
