library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicAlignments)
library(panelcn.mops)

##Sets the bed file for which the intervals for copy number variants are
##calculated
if (file.exists("hg19.bed")){
  export.bed(TxDb.Hsapiens.UCSC.hg19.knownGene, "hg19.bed")
}
bed <- "hg19.bed"
countWindows <- getWindows(bed)

##Labels the folder that contains the sample folders under the path parameter
control_samples <- list.dirs(path = "./healthy", 
                             full.names = TRUE, recursive = FALSE)
experiment_samples <- list.dirs(path = "./patient", 
                                full.names = TRUE, recursive = FALSE)

##Gets the BAM files from sample directories
control_bamfl <- c()
experiment_bamfl <- c()
for (i in control_samples) {
  bam <- list_files_with_exts(dir = i, exts = "bam")
  control_bamfl <- append(control_bamfl, bam)
}

for (i in experiment_samples) {
  bam <- list_files_with_exts(dir = i, exts = "bam")
  experiment_bamfl <- append(experiment_bamfl, bam)
}

##Counts copy number variation of regions found in the bed file
##The results are stored in a dataframe list with each index corresponding to
##each sample
##Control samples are excluded
index <- length(colnames(elementMetadata(experiment)))
control <- countBamListInGRanges(countWindows = countWindows,
                                 bam.files = control_bamfl, 
                                 read.width = 150)
experiment <- countBamListInGRanges(countWindows = countWindows,
                                    bam.files = experiment_bamfl, 
                                    read.width = 150)
elementMetadata(experiment) <- cbind(elementMetadata(experiment),
                                     elementMetadata(control))
resultlist <- runPanelcnMops(experiment, 1:index,
                             countWindows = countWindows)
sampleNames <- colnames(elementMetadata(experiment))
resulttable <- createResultTable(resultlist = resultlist, XandCB = experiment,
                                 countWindows = countWindows,
                                 sampleNames = sampleNames)

##Creates a CSV file containing the copy number counts of all the samples
##CNVRanger has NE_id instead of sample_ID which is more generic
CNV_calls <- data.frame()
for (data in resulttable) {
  data$Sample <- gsub("\\..*","", data$Sample)
  data$Sample <- gsub("_quality_checked_.*","", data$Sample)
  data$CN <- gsub("CN", "", data$CN)
  calls <- data.frame(chr = paste0("chr", data$Chr),
                      start = data$Start,
                      end = data$End,
                      sample_ID = data$Sample,
                      state = data$CN)
  CNV_calls <- merge(CNV_calls, calls, all.x = TRUE, all.y = TRUE)
  CNV_calls <- merge(CNV_calls, calls, all.x = TRUE, all.y = TRUE)
}

write.csv(CNV_calls, "Copy Number Calls.csv", row.names = FALSE)
