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