#' Call single nucleotide polymorphism variants
#' 
#' This function calls SNP variants from BAM files. The results are formatted 
#' into a VCF file and the ID column is populated with dbSNP IDs.
#' 
#' @param bam A list of paths to BAM files
#' @param threads The number of cores to use.
#' @param outputdir The output directory for the VCF file.
#' @return A list of file paths to the VCF files.
#' @export
callsnp <- function(bam,
                    threads = 4L,
                    outputdir = getwd()) {
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
  ##Sets the reference genome that corresponds to the species chosen by the user
  switch(species,
         "1"={
           ##Homo sapiens hg19
           organism <- BSgenome.Hsapiens.UCSC.hg19
           
           ##Selects hg19 as the reference genome
           ##If reference doesn't exist within package directory, create one
           if (dir.exists(paste0(find.package("ExpVar"), "/hg19"))) {
             refgen <- GmapGenome(organism,
                                  directory = find.package("ExpVar"))
           } else {
             print("Reference genome not found. Creating reference. This might take a while...")
             refgen <- GmapGenome(organism,
                                  directory = find.package("ExpVar"),
                                  create = TRUE)
           }
           
           for (i in bam) {
             if (isFALSE(file.exists(paste0(i, ".bai")))) {
               print("Indexing bam...")
               setwd(dirname(i))
               sortBam(i, i)
               indexBam(i)
             }  
             
             setwd(wd)
             bamfl <- BamFile(i)
             vcflist <- c()
             bpp <- BiocParallel::MulticoreParam(threads)
             chromosomes <- seqlevels(refgen)
             tally.Param <- TallyVariantsParam(refgen, high_base_quality = 23L)
             print("Tallying BAM file...")
             tallies <- tallyVariants(bamfl, tally.Param, BPPARAM = bpp)
             gc()
             calling.filters <- VariantCallingFilters()
             post.filters <- VariantPostFilters()
             print("Calling variants...")
             snp <- callVariants(tallies, calling.filters, post.filters)
             sampleNames(snp) <- file_path_sans_ext(basename(i))
             mcols(snp) <- NULL
             print("Formatting variant information as VCF...")
             vcf <- asVCF(snp)
             
             
             vcflist <- c()
             writeVcf(vcf, paste0(file_path_sans_ext(i), ".vcf"),
                      index = TRUE)
             vcf <- readVcf(paste0(file_path_sans_ext(i), ".vcf.bgz"))
             
             ##Due to constraints in memory, rsIDs are obtained on a chromosome by chromosome
             ##basis
             ##The resulting data frame is sorted as per the vcf and the rsIDs are injected
             ##into the vcf before writing it as a file again
             print("Loading dbSNP information...")
             all_snps <- SNPlocs.Hsapiens.dbSNP144.GRCh37
             vcfID <- data.frame()
             mainChromosomes <- paste0("chr", seqlevels(all_snps))
             
             print("Finding rsIDs...")
             for (x in mainChromosomes) {
               vcf_chrom <- vcf[grepl(names(vcf),
                                      pattern = paste0(x, ":"))]
               rd_chr <- rowRanges(vcf_chrom)
               tar_chr <- as.vector(seqnames(rd_chr)@values)
               if (length(tar_chr) > 0) {
                 tar_chr <- gsub("chr", "", tar_chr)
                 tar_chr <- seqlevels(all_snps)[seqlevels(all_snps) == tar_chr]
                 my_snps <- snpsBySeqname(all_snps, c(tar_chr))
                 snp_ID <- data.frame(posIDX = paste0("chr",
                                                      seqnames(my_snps), 
                                                      ":", 
                                                      pos(my_snps)), 
                                      rsID = my_snps$RefSNP_id, 
                                      stringsAsFactors = FALSE)
                 matV1 <- data.frame(Variant = names(rd_chr), stringsAsFactors = FALSE)
                 matV1$chromosome <- gsub("(.*):(.*)_(.*)/(.*)", "\\1", matV1$Variant)
                 matV1$start <- gsub("(.*):(.*)_(.*)/(.*)", "\\2", matV1$Variant)
                 matV1$end <- gsub("(.*):(.*)_(.*)/(.*)", "\\2", matV1$Variant)
                 matV1$ref_allele <- gsub("(.*):(.*)_(.*)/(.*)", "\\3", matV1$Variant)
                 matV1$alt_allele <- gsub("(.*):(.*)_(.*)/(.*)", "\\4", matV1$Variant)
                 matV1$posIDX <- gsub("(.*)_(.*)", "\\1", matV1$Variant)
                 matS <- merge(matV1,
                               snp_ID[!duplicated(snp_ID[, "posIDX"]),],
                               all.x=TRUE,
                               by="posIDX")
                 matS <- dplyr::select(matS,-posIDX)
                 vcfID <- merge(vcfID, matS, all.x = TRUE, all.y = TRUE)
                 vcfID <- merge(vcfID, matS, all.x = TRUE, all.y = TRUE)
               } else {
                 print(paste0(x, " not found..."))
               }
               
             }
             
             print("Injecting rsIDs into VCF...")
             vcfID$chr <-  as.integer(gsub("chr", "", vcfID$chromosome))
             vcfID$chr[vcfID$chromosome == "chrX"] <- 23L
             vcfID$chr[vcfID$chromosome == "chrY"] <- 24L
             vcfID$chr[vcfID$chromosome == "chrM"] <- 25L
             vcfID$chr[vcfID$chromosome == "chrMT"] <- 26L
             dd <- vcfID[order(vcfID$chr, as.integer(vcfID$start)),]
             dd$rsID[is.na(dd$rsID)] <- "."
             rsID <- dd$rsID[dd$chromosome != "chrMT"]
             if (length(vcfID$chr[vcfID$chromosome == "chrMT"]) > 0) {
               M <- rep(".", length(names(vcf[seqnames(vcf) == "chrM",])))
               rsID <- append(rsID, M)
               MT <- dd$rsID[dd$chromosome == "chrMT"]
               rsID <- append(rsID, MT)
             }
             fill <- rep(".", length(names(vcf)) - length(rsID))
             rsID <- append(rsID, fill)
             names(vcf) <- rsID
             
             print("Writing to VCF file...")
             setwd(outputdir)
             a <- writeVcf(vcf, paste0(file_path_sans_ext(i), ".vcf"), 
                           index = TRUE)
             a <- paste0(outputdir, "/", a)
             vcflist <- append(vcflist, a)
           }
         },
         "2"={
           ##Homo sapiens hg38
           organism <- BSgenome.Hsapiens.UCSC.hg38
           
           ##Selects hg38 as the reference genome
           ##If reference doesn't exist within package directory, create one
           if (dir.exists(paste0(find.package("ExpVar"), "/hg38"))) {
             refgen <- GmapGenome(organism,
                                  directory = find.package("ExpVar"))
           } else {
             print("Reference genome not found. Creating reference. This might take a while...")
             refgen <- GmapGenome(organism,
                                  directory = find.package("ExpVar"),
                                  create = TRUE)
           }
           
           for (i in bam) {
             if (isFALSE(file.exists(paste0(i, ".bai")))) {
               print("Indexing bam...")
               setwd(dirname(i))
               sortBam(i, i)
               indexBam(i)
             }  
             
             setwd(wd)
             bamfl <- BamFile(i)
             vcflist <- c()
             bpp <- BiocParallel::MulticoreParam(threads)
             chromosomes <- seqlevels(refgen)
             tally.Param <- TallyVariantsParam(refgen, high_base_quality = 23L)
             print("Tallying BAM file...")
             tallies <- tallyVariants(bamfl, tally.Param, BPPARAM = bpp)
             gc()
             calling.filters <- VariantCallingFilters()
             post.filters <- VariantPostFilters()
             print("Calling variants...")
             snp <- callVariants(tallies, calling.filters, post.filters)
             sampleNames(snp) <- file_path_sans_ext(basename(i))
             mcols(snp) <- NULL
             print("Formatting variant information as VCF...")
             vcf <- asVCF(snp)
             
             
             vcflist <- c()
             writeVcf(vcf, paste0(file_path_sans_ext(i), ".vcf"),
                      index = TRUE)
             vcf <- readVcf(paste0(file_path_sans_ext(i), ".vcf.bgz"))
             
             ##Due to constraints in memory, rsIDs are obtained on a chromosome by chromosome
             ##basis
             ##The resulting data frame is sorted as per the vcf and the rsIDs are injected
             ##into the vcf before writing it as a file again
             print("Loading dbSNP information...")
             all_snps <- SNPlocs.Hsapiens.dbSNP144.GRCh37
             vcfID <- data.frame()
             mainChromosomes <- paste0("chr", seqlevels(all_snps))
             
             print("Finding rsIDs...")
             for (x in mainChromosomes) {
               vcf_chrom <- vcf[grepl(names(vcf),
                                      pattern = paste0(x, ":"))]
               rd_chr <- rowRanges(vcf_chrom)
               tar_chr <- as.vector(seqnames(rd_chr)@values)
               if (length(tar_chr) > 0) {
                 tar_chr <- gsub("chr", "", tar_chr)
                 tar_chr <- seqlevels(all_snps)[seqlevels(all_snps) == tar_chr]
                 my_snps <- snpsBySeqname(all_snps, c(tar_chr))
                 snp_ID <- data.frame(posIDX = paste0("chr",
                                                      seqnames(my_snps), 
                                                      ":", 
                                                      pos(my_snps)), 
                                      rsID = my_snps$RefSNP_id, 
                                      stringsAsFactors = FALSE)
                 matV1 <- data.frame(Variant = names(rd_chr), stringsAsFactors = FALSE)
                 matV1$chromosome <- gsub("(.*):(.*)_(.*)/(.*)", "\\1", matV1$Variant)
                 matV1$start <- gsub("(.*):(.*)_(.*)/(.*)", "\\2", matV1$Variant)
                 matV1$end <- gsub("(.*):(.*)_(.*)/(.*)", "\\2", matV1$Variant)
                 matV1$ref_allele <- gsub("(.*):(.*)_(.*)/(.*)", "\\3", matV1$Variant)
                 matV1$alt_allele <- gsub("(.*):(.*)_(.*)/(.*)", "\\4", matV1$Variant)
                 matV1$posIDX <- gsub("(.*)_(.*)", "\\1", matV1$Variant)
                 matS <- merge(matV1,
                               snp_ID[!duplicated(snp_ID[, "posIDX"]),],
                               all.x=TRUE,
                               by="posIDX")
                 matS <- dplyr::select(matS,-posIDX)
                 vcfID <- merge(vcfID, matS, all.x = TRUE, all.y = TRUE)
                 vcfID <- merge(vcfID, matS, all.x = TRUE, all.y = TRUE)
               } else {
                 print(paste0(x, " not found..."))
               }
               
             }
             
             print("Injecting rsIDs into VCF...")
             vcfID$chr <-  as.integer(gsub("chr", "", vcfID$chromosome))
             vcfID$chr[vcfID$chromosome == "chrX"] <- 23L
             vcfID$chr[vcfID$chromosome == "chrY"] <- 24L
             vcfID$chr[vcfID$chromosome == "chrM"] <- 25L
             vcfID$chr[vcfID$chromosome == "chrMT"] <- 26L
             dd <- vcfID[order(vcfID$chr, as.integer(vcfID$start)),]
             dd$rsID[is.na(dd$rsID)] <- "."
             rsID <- dd$rsID[dd$chromosome != "chrMT"]
             if (length(vcfID$chr[vcfID$chromosome == "chrMT"]) > 0) {
               M <- rep(".", length(names(vcf[seqnames(vcf) == "chrM",])))
               rsID <- append(rsID, M)
               MT <- dd$rsID[dd$chromosome == "chrMT"]
               rsID <- append(rsID, MT)
             }
             fill <- rep(".", length(names(vcf)) - length(rsID))
             rsID <- append(rsID, fill)
             names(vcf) <- rsID
             
             print("Writing to VCF file...")
             setwd(outputdir)
             a <- writeVcf(vcf, paste0(file_path_sans_ext(i), ".vcf"), 
                           index = TRUE)
             a <- paste0(outputdir, "/", a)
             vcflist <- append(vcflist, a)
           }
         },
         "3"={
           ##Mus musculus mm10
           organism <- BSgenome.Mmusculus.UCSC.mm10
           
           ##Selects mm10 as the reference genome
           ##If reference doesn't exist within package directory, create one
           if (dir.exists(paste0(find.package("ExpVar"), "/mm10"))) {
             refgen <- GmapGenome(organism,
                                  directory = find.package("ExpVar"))
           } else {
             print("Reference genome not found. Creating reference. This might take a while...")
             refgen <- GmapGenome(organism,
                                  directory = find.package("ExpVar"),
                                  create = TRUE)
           }
           
           for (i in bam) {
             if (isFALSE(file.exists(paste0(i, ".bai")))) {
               print("Indexing bam...")
               setwd(dirname(i))
               sortBam(i, i)
               indexBam(i)
             }  
             
             setwd(wd)
             bamfl <- BamFile(i)
             vcflist <- c()
             bpp <- BiocParallel::MulticoreParam(threads)
             chromosomes <- seqlevels(refgen)
             tally.Param <- TallyVariantsParam(refgen, high_base_quality = 23L)
             print("Tallying BAM file...")
             tallies <- tallyVariants(bamfl, tally.Param, BPPARAM = bpp)
             gc()
             calling.filters <- VariantCallingFilters()
             post.filters <- VariantPostFilters()
             print("Calling variants...")
             snp <- callVariants(tallies, calling.filters, post.filters)
             sampleNames(snp) <- file_path_sans_ext(basename(i))
             mcols(snp) <- NULL
             print("Formatting variant information as VCF...")
             vcf <- asVCF(snp)
             
             
             vcflist <- c()
             print("Writing to VCF file...")
             setwd(outputdir)
             a <- writeVcf(vcf, paste0(file_path_sans_ext(i), ".vcf"), 
                           index = TRUE)
             a <- paste0(outputdir, "/", a)
             vcflist <- append(vcflist, a)
         }},
         "4"={
           ##Arabidopsis thaliana TAIR9
           organism <- BSgenome.Athaliana.TAIR.TAIR9
           
           ##Selects hg19 as the reference genome
           ##If reference doesn't exist within package directory, create one
           if (dir.exists(paste0(find.package("ExpVar"), "/TAIR9"))) {
             refgen <- GmapGenome(organism,
                                  directory = find.package("ExpVar"))
           } else {
             print("Reference genome not found. Creating reference. This might take a while...")
             refgen <- GmapGenome(organism,
                                  directory = find.package("ExpVar"),
                                  create = TRUE)
           }
           
           for (i in bam) {
             if (isFALSE(file.exists(paste0(i, ".bai")))) {
               print("Indexing bam...")
               setwd(dirname(i))
               sortBam(i, i)
               indexBam(i)
             }  
             
             setwd(wd)
             bamfl <- BamFile(i)
             vcflist <- c()
             bpp <- BiocParallel::MulticoreParam(threads)
             chromosomes <- seqlevels(refgen)
             tally.Param <- TallyVariantsParam(refgen, high_base_quality = 23L)
             print("Tallying BAM file...")
             tallies <- tallyVariants(bamfl, tally.Param, BPPARAM = bpp)
             gc()
             calling.filters <- VariantCallingFilters()
             post.filters <- VariantPostFilters()
             print("Calling variants...")
             snp <- callVariants(tallies, calling.filters, post.filters)
             sampleNames(snp) <- file_path_sans_ext(basename(i))
             mcols(snp) <- NULL
             print("Formatting variant information as VCF...")
             vcf <- asVCF(snp)
             
             
             vcflist <- c()
             print("Writing to VCF file...")
             setwd(outputdir)
             a <- writeVcf(vcf, paste0(file_path_sans_ext(i), ".vcf"), 
                           index = TRUE)
             a <- paste0(outputdir, "/", a)
             vcflist <- append(vcflist, a)
         }},
         "5"={
           ##Drosophilia melanogaster dm6
           organism <- BSgenome.Dmelanogaster.UCSC.dm6
           
           ##Selects dm6 as the reference genome
           ##If reference doesn't exist within package directory, create one
           if (dir.exists(paste0(find.package("ExpVar"), "/dm6"))) {
             refgen <- GmapGenome(organism,
                                  directory = find.package("ExpVar"))
           } else {
             print("Reference genome not found. Creating reference. This might take a while...")
             refgen <- GmapGenome(organism,
                                  directory = find.package("ExpVar"),
                                  create = TRUE)
           }
           
           for (i in bam) {
             if (isFALSE(file.exists(paste0(i, ".bai")))) {
               print("Indexing bam...")
               setwd(dirname(i))
               sortBam(i, i)
               indexBam(i)
             }  
             
             setwd(wd)
             bamfl <- BamFile(i)
             vcflist <- c()
             bpp <- BiocParallel::MulticoreParam(threads)
             chromosomes <- seqlevels(refgen)
             tally.Param <- TallyVariantsParam(refgen, high_base_quality = 23L)
             print("Tallying BAM file...")
             tallies <- tallyVariants(bamfl, tally.Param, BPPARAM = bpp)
             gc()
             calling.filters <- VariantCallingFilters()
             post.filters <- VariantPostFilters()
             print("Calling variants...")
             snp <- callVariants(tallies, calling.filters, post.filters)
             sampleNames(snp) <- file_path_sans_ext(basename(i))
             mcols(snp) <- NULL
             print("Formatting variant information as VCF...")
             vcf <- asVCF(snp)
             
             
             vcflist <- c()
             print("Writing to VCF file...")
             setwd(outputdir)
             a <- writeVcf(vcf, paste0(file_path_sans_ext(i), ".vcf"), 
                           index = TRUE)
             a <- paste0(outputdir, "/", a)
             vcflist <- append(vcflist, a)
         }},
         "6"={
           ##Danio rerio danRer11
           organism <- BSgenome.Drerio.UCSC.danRer11
           
           ##Selects danRer11 as the reference genome
           ##If reference doesn't exist within package directory, create one
           if (dir.exists(paste0(find.package("ExpVar"), "/danRer11"))) {
             refgen <- GmapGenome(organism,
                                  directory = find.package("ExpVar"))
           } else {
             print("Reference genome not found. Creating reference. This might take a while...")
             refgen <- GmapGenome(organism,
                                  directory = find.package("ExpVar"),
                                  create = TRUE)
           }
           
           for (i in bam) {
             if (isFALSE(file.exists(paste0(i, ".bai")))) {
               print("Indexing bam...")
               setwd(dirname(i))
               sortBam(i, i)
               indexBam(i)
             }  
             
             setwd(wd)
             bamfl <- BamFile(i)
             vcflist <- c()
             bpp <- BiocParallel::MulticoreParam(threads)
             chromosomes <- seqlevels(refgen)
             tally.Param <- TallyVariantsParam(refgen, high_base_quality = 23L)
             print("Tallying BAM file...")
             tallies <- tallyVariants(bamfl, tally.Param, BPPARAM = bpp)
             gc()
             calling.filters <- VariantCallingFilters()
             post.filters <- VariantPostFilters()
             print("Calling variants...")
             snp <- callVariants(tallies, calling.filters, post.filters)
             sampleNames(snp) <- file_path_sans_ext(basename(i))
             mcols(snp) <- NULL
             print("Formatting variant information as VCF...")
             vcf <- asVCF(snp)
             
             
             vcflist <- c()
             print("Writing to VCF file...")
             setwd(outputdir)
             a <- writeVcf(vcf, paste0(file_path_sans_ext(i), ".vcf"), 
                           index = TRUE)
             a <- paste0(outputdir, "/", a)
             vcflist <- append(vcflist, a)
         }},
         "7"={
           ##Rattus norvegicus rn5
           organism <- BSgenome.Rnorvegicus.UCSC.rn5
           
           ##Selects danRer11 as the reference genome
           ##If reference doesn't exist within package directory, create one
           if (dir.exists(paste0(find.package("ExpVar"), "/rn5"))) {
             refgen <- GmapGenome(organism,
                                  directory = find.package("ExpVar"))
           } else {
             print("Reference genome not found. Creating reference. This might take a while...")
             refgen <- GmapGenome(organism,
                                  directory = find.package("ExpVar"),
                                  create = TRUE)
           }
           
           for (i in bam) {
             if (isFALSE(file.exists(paste0(i, ".bai")))) {
               print("Indexing bam...")
               setwd(dirname(i))
               sortBam(i, i)
               indexBam(i)
             }  
             
             setwd(wd)
             bamfl <- BamFile(i)
             vcflist <- c()
             bpp <- BiocParallel::MulticoreParam(threads)
             chromosomes <- seqlevels(refgen)
             tally.Param <- TallyVariantsParam(refgen, high_base_quality = 23L)
             print("Tallying BAM file...")
             tallies <- tallyVariants(bamfl, tally.Param, BPPARAM = bpp)
             gc()
             calling.filters <- VariantCallingFilters()
             post.filters <- VariantPostFilters()
             print("Calling variants...")
             snp <- callVariants(tallies, calling.filters, post.filters)
             sampleNames(snp) <- file_path_sans_ext(basename(i))
             mcols(snp) <- NULL
             print("Formatting variant information as VCF...")
             vcf <- asVCF(snp)
             
             
             vcflist <- c()
             print("Writing to VCF file...")
             setwd(outputdir)
             a <- writeVcf(vcf, paste0(file_path_sans_ext(i), ".vcf"), 
                           index = TRUE)
             a <- paste0(outputdir, "/", a)
             vcflist <- append(vcflist, a)
         }},
         "8"={
           ##Saccharomyces cerevisiae sacCer3
           organism <- BSgenome.Scerevisiae.UCSC.sacCer3
           
           ##Selects sacCer3 as the reference genome
           ##If reference doesn't exist within package directory, create one
           if (dir.exists(paste0(find.package("ExpVar"), "/sacCer3"))) {
             refgen <- GmapGenome(organism,
                                  directory = find.package("ExpVar"))
           } else {
             print("Reference genome not found. Creating reference. This might take a while...")
             refgen <- GmapGenome(organism,
                                  directory = find.package("ExpVar"),
                                  create = TRUE)
           }
           
           for (i in bam) {
             if (isFALSE(file.exists(paste0(i, ".bai")))) {
               print("Indexing bam...")
               setwd(dirname(i))
               sortBam(i, i)
               indexBam(i)
             }  
             
             setwd(wd)
             bamfl <- BamFile(i)
             vcflist <- c()
             bpp <- BiocParallel::MulticoreParam(threads)
             chromosomes <- seqlevels(refgen)
             tally.Param <- TallyVariantsParam(refgen, high_base_quality = 23L)
             print("Tallying BAM file...")
             tallies <- tallyVariants(bamfl, tally.Param, BPPARAM = bpp)
             gc()
             calling.filters <- VariantCallingFilters()
             post.filters <- VariantPostFilters()
             print("Calling variants...")
             snp <- callVariants(tallies, calling.filters, post.filters)
             sampleNames(snp) <- file_path_sans_ext(basename(i))
             mcols(snp) <- NULL
             print("Formatting variant information as VCF...")
             vcf <- asVCF(snp)
             
             
             vcflist <- c()
             print("Writing to VCF file...")
             setwd(outputdir)
             a <- writeVcf(vcf, paste0(file_path_sans_ext(i), ".vcf"), 
                           index = TRUE)
             a <- paste0(outputdir, "/", a)
             vcflist <- append(vcflist, a)
         }},
         "9"={
           ##Caenorhabditis elagans
           organism <- BSgenome.Celegans.UCSC.ce11
           
           ##Selects ce11 as the reference genome
           ##If reference doesn't exist within package directory, create one
           if (dir.exists(paste0(find.package("ExpVar"), "/ce11"))) {
             refgen <- GmapGenome(organism,
                                  directory = find.package("ExpVar"))
           } else {
             print("Reference genome not found. Creating reference. This might take a while...")
             refgen <- GmapGenome(organism,
                                  directory = find.package("ExpVar"),
                                  create = TRUE)
           }
           
           for (i in bam) {
             if (isFALSE(file.exists(paste0(i, ".bai")))) {
               print("Indexing bam...")
               setwd(dirname(i))
               sortBam(i, i)
               indexBam(i)
             }  
             
             setwd(wd)
             bamfl <- BamFile(i)
             vcflist <- c()
             bpp <- BiocParallel::MulticoreParam(threads)
             chromosomes <- seqlevels(refgen)
             tally.Param <- TallyVariantsParam(refgen, high_base_quality = 23L)
             print("Tallying BAM file...")
             tallies <- tallyVariants(bamfl, tally.Param, BPPARAM = bpp)
             gc()
             calling.filters <- VariantCallingFilters()
             post.filters <- VariantPostFilters()
             print("Calling variants...")
             snp <- callVariants(tallies, calling.filters, post.filters)
             sampleNames(snp) <- file_path_sans_ext(basename(i))
             mcols(snp) <- NULL
             print("Formatting variant information as VCF...")
             vcf <- asVCF(snp)
             
             
             vcflist <- c()
             print("Writing to VCF file...")
             setwd(outputdir)
             a <- writeVcf(vcf, paste0(file_path_sans_ext(i), ".vcf"), 
                           index = TRUE)
             a <- paste0(outputdir, "/", a)
             vcflist <- append(vcflist, a)
           }
         }
  )
             setwd(wd)
             return(vcflist)
}


