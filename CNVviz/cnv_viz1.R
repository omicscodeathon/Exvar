#Reading and accessing CNV data

call.file <- file.path(input$csv_file ) 
calls <- read.csv(call.file, as.is=TRUE)
nrow(calls)
head(calls)
length(unique(calls[,"NE_id"]))

#Representation as a GRangesList  
grl <- makeGRangesListFromDataFrame(calls, 
                                    split.field = "NE_id", keep.extra.columns = TRUE)
# sorting of the calls according to their genomic coordinates.
grl <- sort(grl)
#-------------------------------------
#Representation as a RaggedExperiment
ra <- RaggedExperiment(grl)
assay(ra[1:2,1:2])
weight <- rnorm(ncol(ra), mean=1100, sd=100)
fcr <- rnorm(ncol(ra), mean=6, sd=1.5)
colData(ra)$weight <- round(weight, digits=2)
colData(ra)$fcr <- round(fcr, digits=2)
colData(ra)
#-------------------------------------
#Summarizing individual CNV calls across a population 

##Trimming low-density areas
cnvrs <- populationRanges(grl, density = 0.1)
cnvrs <- populationRanges(grl, density = 0.1, est.recur = TRUE)
#and filter for recurrent CNVs that exceed a significance threshold ....
pval = 0.05  # = default = demo
pval <- input$pvalue  # = input
subset(cnvrs, pvalue < pval)    

# illustrate the landscape of recurrent CNV regions 
chr = "chr1" #default = demo
chr <- input$chr
 
output$recurrent_plot <- renderPlot({
  plotRecurrentRegions(cnvrs, genome = "hg38", chr = chr) 
})

#-------------------------------------
#Overlap analysis of CNVs with functional genomic regions  
ah <- AnnotationHub::AnnotationHub()  # request approval for directory creation
ahDb <- AnnotationHub::query(ah, pattern = c("Homo sapiens", "EnsDb"))

# retrieve gene coordinates 
ahEdb <- ahDb[["AH53211"]]

bt.genes <- ensembldb::genes(ahEdb)
GenomeInfoDb::seqlevelsStyle(bt.genes) <- "UCSC"
bt.genes

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


#illustrate 
output$overlap_plot <- renderPlot({
  cnvOncoPrint(grl, cgenes)
})
#-------------------------------------
#Overlap permutation test 

res <- suppressWarnings(overlapPermTest(A=sel.cnvrs, B=sel.genes, ntimes=100, 
                                        genome="hg19", mask=NA, per.chromosome=TRUE, count.once=TRUE))

summary(res[[1]]$permuted)

#plot
output$permutation_plot <- renderPlot({
  plot(res)
})













