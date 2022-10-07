# data

data.dir <- system.file("extdata", package="CNVRanger")
call.file <- file.path(data.dir, "Silva16_PONE_CNV_calls.csv")
calls <- read.csv(call.file, as.is=TRUE)
nrow(calls)
length(unique(calls[,"NE_id"]))
#Representation as a GRangesList
grl <- makeGRangesListFromDataFrame(calls, 
                                    split.field = "NE_id", keep.extra.columns = TRUE)
grl <- sort(grl)
grl

#Representation as a RaggedExperiment
ra <- RaggedExperiment(grl)
ra

assay(ra[1:5,1:5])

weight <- rnorm(ncol(ra), mean=1100, sd=100)
fcr <- rnorm(ncol(ra), mean=6, sd=1.5)
colData(ra)$weight <- round(weight, digits=2)
colData(ra)$fcr <- round(fcr, digits=2)
colData(ra)

#Trimming low-density areas
cnvrs <- populationRanges(grl, density = 0.1)
cnvrs

#Identifying recurrent regions
cnvrs <- populationRanges(grl, density = 0.1, est.recur = TRUE)
cnvrs
pval = 0.05  # = default = demo
pval <- input$pvalue  # = input
subset(cnvrs, pvalue < pval)  
#
chr = "chr1" #default = demo
chr <- input$chr
# plot1
output$recurrent_plot <- renderPlot({
  plotRecurrentRegions(cnvrs, genome = "bosTau6", chr = chr)
})
#####-----
#Overlap analysis of CNVs with functional genomic regions

ah <- AnnotationHub::AnnotationHub()
ahDb <- AnnotationHub::query(ah, pattern = c("Bos taurus", "EnsDb"))
ahDb
ahEdb <- ahDb[["AH60948"]]
bt.genes <- ensembldb::genes(ahEdb)
GenomeInfoDb::seqlevelsStyle(bt.genes) <- "UCSC"
bt.genes

# filter sequence = iput

sel.genes <- subset(bt.genes, seqnames %in% paste0("chr", 1:2))
sel.genes <- subset(sel.genes, gene_biotype == "protein_coding") 
sel.cnvrs <- subset(cnvrs, seqnames %in% paste0("chr", 1:2))

#Finding and illustrating overlaps
olaps <- GenomicRanges::findOverlaps(sel.genes, sel.cnvrs, ignore.strand=TRUE)
qh <- S4Vectors::queryHits(olaps)
sh <- S4Vectors::subjectHits(olaps)
cgenes <- sel.genes[qh]
cgenes$type <- sel.cnvrs$type[sh]
subset(cgenes, select = "type")

# plot 2
output$overlap_plot <- renderPlot({
  cnvOncoPrint(grl, cgenes)
})


####---------------------
#Overlap permutation test

res <- suppressWarnings(overlapPermTest(A=sel.cnvrs, B=sel.genes, ntimes=100, 
                                        genome="bosTau6", mask=NA, per.chromosome=TRUE, count.once=TRUE))
res
summary(res[[1]]$permuted)

#plot3
output$permutation_plot <- renderPlot({
  plot(res)
})



