#' Call single nucleotide polymorphism variants
#' 
#' This function calls SNP variants from BAM files. The results are formatted 
#' into a VCF file and the ID column is populated with dbSNP IDs (Homo sapiens
#' only).
#' 
#' @param bam A list of paths to BAM files
#' @param threads The number of cores to use.
#' @param outputdir The output directory for the VCF file.
#' @return A list of file paths to the VCF files.
#' @export
callsnp <- function(bam,
                    threads = 4L,
                    outputdir = getwd()) {
    if(Sys.info()[['sysname']] != "Linux"){
    message("This function is only available on Linux.")
    stop()
    }
  cat(paste0("These are the species currently supported by Exvar: \n",
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
  cat("Enable single-chromosome splitting? \nNote: enabling splitting reduces the RAM requirements, but may take longer.")
  mode <- readline("Type [y/n] for [yes/no]: ")
  wd <- getwd()
  
  library(gmapR)
  library(Rsamtools)
  library(VariantTools)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(VariantAnnotation)
  library(BiocParallel)
  
  if(mode == "y") {
    ##Sets the reference genome that corresponds to the species chosen by the user
    switch(species,
           "1"={
             ##Homo sapiens hg19
             library(BSgenome.Hsapiens.UCSC.hg19)
             library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
             organism <- BSgenome.Hsapiens.UCSC.hg19
             
             ##Selects hg19 as the reference genome
             ##If reference doesn't exist within package directory, create one
             if (dir.exists(paste0(find.package("Exvar"), "/hg19"))) {
               refgen <- GmapGenome(organism,
                                    directory = find.package("Exvar"))
             } else {
               print("Reference genome not found. Creating reference. This might take a while...")
               refgen <- GmapGenome(organism,
                                    directory = find.package("Exvar"),
                                    create = TRUE)
             }
             output <- c()
             for (i in bam) {
               print(paste0("Analysing ", basename(i)))
               if (isFALSE(file.exists(paste0(i, ".bai")))) {
                 print("Indexing bam...")
                 setwd(dirname(i))
                 sortBam(i, file_path_sans_ext(i))
                 indexBam(i)
               }  
               setwd(wd)
               bamfl <- BamFile(i)
               chromosomes <- seqlevels(bamfl)
               print("Splitting BAM into chromosomes...")
               tempfolder <- paste0(dirname(i), "/temp")
               dir.create(tempfolder)
               for (x in chromosomes) {
                 print(paste0("Creating ", 
                              basename(file_path_sans_ext(i)), "_", x, ".bam ",
                              "in temp folder."))
                 system(paste0("samtools view -b ", 
                               i,
                               " ", x, " > ",
                               tempfolder, "/", 
                               basename(file_path_sans_ext(i)), "_", x, ".bam"))
                 setwd(tempfolder)
                 sortBam(paste0(basename(file_path_sans_ext(i)), "_", x, ".bam"),
                         paste0(basename(file_path_sans_ext(i)), "_", x))
                 indexBam(paste0(basename(file_path_sans_ext(i)), "_", x, ".bam"))
                 setwd(wd)
               }
               chrbam <- list_files_with_exts(tempfolder, "bam")
               chrbam <- sort(chrbam)
               
               
               for (bam in chrbam) {
                 print(paste0("Identifying variants in ", basename(bam)))
                 bamfl <- BamFile(bam)
                 bpp <- BiocParallel::MulticoreParam(threads)
                 scan <- scanBam(bam)
                 if (length(scan[[1]]$rname) < 1) next
                 chromosome <- as.character(scan[[1]]$rname[[1]])
                 rm(scan)
                 gr <- GRanges(seqinfo(refgen))
                 gr <- gr[seqnames(gr) != chromosome]
                 tally.Param <- TallyVariantsParam(refgen, high_base_quality = 23L,
                                                   mask = gr)
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
                 
                 
                 
                 vcfpath <- c()
                 if (length(snp) > 0) {
                   vcf <- asVCF(snp)
                   writeVcf(vcf, paste0(file_path_sans_ext(bam), ".vcf"),
                            index = FALSE)
                   vcf <- readVcf(paste0(file_path_sans_ext(bam), ".vcf"))
                   
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
                     }
                   }
                   
                   print("Writing to VCF file...")
                   setwd(tempfolder)
                   print(getwd())
                   print(paste0(file_path_sans_ext(basename(bam)), ".vcf"))
                   a <- writeVcf(vcf, paste0(file_path_sans_ext(basename(bam)), ".vcf"), 
                                 index = FALSE)
                   gc()
                   setwd(wd)
                 } else {
                   print(paste0("No variants found in ", (basename(bam))))
                 }
               }
               
               print("Merging VCF files...")
               vcflist <- list_files_with_exts(tempfolder, "vcf")
               vcfpath <- normalizePath(vcflist)
               match <- VariantAnnotation::readVcf(vcfpath[1])
               print(paste0(gsub(".vcf", "", basename(vcfpath[1])), " OKAY"))
               for (item in vcfpath[-1]) {
                 b <- VariantAnnotation::readVcf(item)
                 match <- BiocGenerics::rbind(match, b)
                 print(paste0(gsub(".vcf", "", basename(item)), " OKAY"))
               }
               rsIDs <- rownames(match)
               rownames(match) <- gsub("chr?(.*)", ".", rsIDs)
               writeVcf(match, paste0(outputdir, "/",
                                      basename(file_path_sans_ext(i)),
                                      "_SNP", ".vcf"), 
                        index = TRUE)
               file.remove(list.files(normalizePath(tempfolder), full.names = TRUE))
               print(paste0("Created ", file_path_sans_ext(basename(i)), "_SNP", ".vcf.bgz"))
               gc()
               setwd(wd)
               append(output, paste0(outputdir, "/", 
                                     file_path_sans_ext(basename(i)),
                                     "_SNP", ".vcf.bgz"))
             }
             file.remove(normalizePath(tempfolder))
           },
           "2"={
             ##Homo sapiens hg38
             library(BSgenome.Hsapiens.UCSC.hg38)
             library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
             organism <- BSgenome.Hsapiens.UCSC.hg38
             
             ##Selects hg38 as the reference genome
             ##If reference doesn't exist within package directory, create one
             if (dir.exists(paste0(find.package("Exvar"), "/hg38"))) {
               refgen <- GmapGenome(organism,
                                    directory = find.package("Exvar"))
             } else {
               print("Reference genome not found. Creating reference. This might take a while...")
               refgen <- GmapGenome(organism,
                                    directory = find.package("Exvar"),
                                    create = TRUE)
             }
             output <- c()
             for (i in bam) {
               print(paste0("Analysing ", basename(i)))
               if (isFALSE(file.exists(paste0(i, ".bai")))) {
                 print("Indexing bam...")
                 setwd(dirname(i))
                 sortBam(i, file_path_sans_ext(i))
                 indexBam(i)
               }  
               setwd(wd)
               bamfl <- BamFile(i)
               chromosomes <- seqlevels(bamfl)
               print("Splitting BAM into chromosomes...")
               tempfolder <- paste0(dirname(i), "/temp")
               dir.create(tempfolder)
               for (x in chromosomes) {
                 print(paste0("Creating ", 
                              basename(file_path_sans_ext(i)), "_", x, ".bam ",
                              "in temp folder."))
                 system(paste0("samtools view -b ", 
                               i,
                               " ", x, " > ",
                               tempfolder, "/", 
                               basename(file_path_sans_ext(i)), "_", x, ".bam"))
                 setwd(tempfolder)
                 sortBam(paste0(basename(file_path_sans_ext(i)), "_", x, ".bam"),
                         paste0(basename(file_path_sans_ext(i)), "_", x))
                 indexBam(paste0(basename(file_path_sans_ext(i)), "_", x, ".bam"))
                 setwd(wd)
               }
               chrbam <- list_files_with_exts(tempfolder, "bam")
               chrbam <- sort(chrbam)
               
               
               for (bam in chrbam) {
                 print(paste0("Identifying variants in ", basename(bam)))
                 bamfl <- BamFile(bam)
                 bpp <- BiocParallel::MulticoreParam(threads)
                 scan <- scanBam(bam)
                 if (length(scan[[1]]$rname) < 1) next
                 chromosome <- as.character(scan[[1]]$rname[[1]])
                 rm(scan)
                 gr <- GRanges(seqinfo(refgen))
                 gr <- gr[seqnames(gr) != chromosome]
                 tally.Param <- TallyVariantsParam(refgen, high_base_quality = 23L,
                                                   mask = gr)
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
                 
                 
                 
                 vcfpath <- c()
                 if (length(snp) > 0) {
                   vcf <- asVCF(snp)
                   writeVcf(vcf, paste0(file_path_sans_ext(bam), ".vcf"),
                            index = FALSE)
                   vcf <- readVcf(paste0(file_path_sans_ext(bam), ".vcf"))
                   
                   ##Due to constraints in memory, rsIDs are obtained on a chromosome by chromosome
                   ##basis
                   ##The resulting data frame is sorted as per the vcf and the rsIDs are injected
                   ##into the vcf before writing it as a file again
                   print("Loading dbSNP information...")
                   all_snps <- SNPlocs.Hsapiens.dbSNP151.GRCh38
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
                     }
                   }
                   
                   print("Writing to VCF file...")
                   setwd(tempfolder)
                   print(getwd())
                   print(paste0(file_path_sans_ext(basename(bam)), ".vcf"))
                   a <- writeVcf(vcf, paste0(file_path_sans_ext(basename(bam)), ".vcf"), 
                                 index = FALSE)
                   gc()
                   setwd(wd)
                 } else {
                   print(paste0("No variants found in ", (basename(bam))))
                 }
               }
               
               print("Merging VCF files...")
               vcflist <- list_files_with_exts(tempfolder, "vcf")
               vcfpath <- normalizePath(vcflist)
               match <- VariantAnnotation::readVcf(vcfpath[1])
               print(paste0(gsub(".vcf", "", basename(vcfpath[1])), " OKAY"))
               for (item in vcfpath[-1]) {
                 b <- VariantAnnotation::readVcf(item)
                 match <- BiocGenerics::rbind(match, b)
                 print(paste0(gsub(".vcf", "", basename(item)), " OKAY"))
               }
               rsIDs <- rownames(match)
               rownames(match) <- gsub("chr?(.*)", ".", rsIDs)
               writeVcf(match, paste0(outputdir, "/",
                                      basename(file_path_sans_ext(i)),
                                      "_SNP", ".vcf"), 
                        index = TRUE)
               file.remove(list.files(normalizePath(tempfolder), full.names = TRUE))
               print(paste0("Created ", file_path_sans_ext(basename(i)), "_SNP", ".vcf.bgz"))
               gc()
               setwd(wd)
               append(output, paste0(outputdir, "/", 
                                     file_path_sans_ext(basename(i)),
                                     "_SNP", ".vcf.bgz"))
             }
             file.remove(normalizePath(tempfolder))
           },
           "3"={
             ##Mus musculus mm10
             library(BSgenome.Mmusculus.UCSC.mm10)
             organism <- BSgenome.Mmusculus.UCSC.mm10
             
             ##Selects mm10 as the reference genome
             ##If reference doesn't exist within package directory, create one
             if (dir.exists(paste0(find.package("Exvar"), "/mm10"))) {
               refgen <- GmapGenome(organism,
                                    directory = find.package("Exvar"))
             } else {
               print("Reference genome not found. Creating reference. This might take a while...")
               refgen <- GmapGenome(organism,
                                    directory = find.package("Exvar"),
                                    create = TRUE)
             }
             output <- c()
             for (i in bam) {
               print(paste0("Analysing ", basename(i)))
               if (isFALSE(file.exists(paste0(i, ".bai")))) {
                 print("Indexing bam...")
                 setwd(dirname(i))
                 sortBam(i, file_path_sans_ext(i))
                 indexBam(i)
               }  
               setwd(wd)
               bamfl <- BamFile(i)
               chromosomes <- seqlevels(bamfl)
               print("Splitting BAM into chromosomes...")
               tempfolder <- paste0(dirname(i), "/temp")
               dir.create(tempfolder)
               for (x in chromosomes) {
                 print(paste0("Creating ", 
                              basename(file_path_sans_ext(i)), "_", x, ".bam ",
                              "in temp folder."))
                 system(paste0("samtools view -b ", 
                               i,
                               " ", x, " > ",
                               tempfolder, "/", 
                               basename(file_path_sans_ext(i)), "_", x, ".bam"))
                 setwd(tempfolder)
                 sortBam(paste0(basename(file_path_sans_ext(i)), "_", x, ".bam"),
                         paste0(basename(file_path_sans_ext(i)), "_", x))
                 indexBam(paste0(basename(file_path_sans_ext(i)), "_", x, ".bam"))
                 setwd(wd)
               }
               chrbam <- list_files_with_exts(tempfolder, "bam")
               chrbam <- sort(chrbam)
               
               
               for (bam in chrbam) {
                 print(paste0("Identifying variants in ", basename(bam)))
                 bamfl <- BamFile(bam)
                 bpp <- BiocParallel::MulticoreParam(threads)
                 scan <- scanBam(bam)
                 if (length(scan[[1]]$rname) < 1) next
                 chromosome <- as.character(scan[[1]]$rname[[1]])
                 rm(scan)
                 gr <- GRanges(seqinfo(refgen))
                 gr <- gr[seqnames(gr) != chromosome]
                 tally.Param <- TallyVariantsParam(refgen, high_base_quality = 23L,
                                                   mask = gr)
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
                 
                 if (length(snp) > 0)  {
                   vcf <- asVCF(snp)
                   print("Writing to VCF file...")
                   setwd(tempfolder)
                   print(getwd())
                   print(paste0(file_path_sans_ext(basename(bam)), ".vcf"))
                   a <- writeVcf(vcf, paste0(file_path_sans_ext(basename(bam)), ".vcf"), 
                                 index = FALSE)
                   gc()
                 } else {
                   print(paste0("No variants found in ", (basename(bam))))
                 }
                 
                 setwd(wd)
               } 
               
               print("Merging VCF files...")
               vcflist <- list_files_with_exts(tempfolder, "vcf")
               vcfpath <- normalizePath(vcflist)
               match <- VariantAnnotation::readVcf(vcfpath[1])
               print(paste0(gsub(".vcf", "", basename(vcfpath[1])), " OKAY"))
               for (item in vcfpath[-1]) {
                 b <- VariantAnnotation::readVcf(item)
                 match <- BiocGenerics::rbind(match, b)
                 print(paste0(gsub(".vcf", "", basename(item)), " OKAY"))
               }
               rsIDs <- rownames(match)
               rownames(match) <- gsub("chr?(.*)", ".", rsIDs)
               writeVcf(match, paste0(outputdir, "/",
                                      basename(file_path_sans_ext(i)),
                                      "_SNP", ".vcf"), 
                        index = TRUE)
               file.remove(list.files(normalizePath(tempfolder), full.names = TRUE))
               print(paste0("Created ", file_path_sans_ext(basename(i)), "_SNP", ".vcf.bgz"))
               gc()
               setwd(wd)
               append(output, paste0(outputdir, "/", 
                                     file_path_sans_ext(basename(i)),
                                     "_SNP", ".vcf.bgz"))
             }
             file.remove(normalizePath(tempfolder))
           },
           "4"={
             ##Arabidopsis thaliana TAIR9
             library(BSgenome.Athaliana.TAIR.TAIR9)
             organism <- BSgenome.Athaliana.TAIR.TAIR9
             
             ##Selects hg19 as the reference genome
             ##If reference doesn't exist within package directory, create one
             if (dir.exists(paste0(find.package("Exvar"), "/TAIR9"))) {
               refgen <- GmapGenome(organism,
                                    directory = find.package("Exvar"))
             } else {
               print("Reference genome not found. Creating reference. This might take a while...")
               refgen <- GmapGenome(organism,
                                    directory = find.package("Exvar"),
                                    create = TRUE)
             }
             output <- c()
             for (i in bam) {
               print(paste0("Analysing ", basename(i)))
               if (isFALSE(file.exists(paste0(i, ".bai")))) {
                 print("Indexing bam...")
                 setwd(dirname(i))
                 sortBam(i, file_path_sans_ext(i))
                 indexBam(i)
               }  
               setwd(wd)
               bamfl <- BamFile(i)
               chromosomes <- seqlevels(bamfl)
               print("Splitting BAM into chromosomes...")
               tempfolder <- paste0(dirname(i), "/temp")
               dir.create(tempfolder)
               for (x in chromosomes) {
                 print(paste0("Creating ", 
                              basename(file_path_sans_ext(i)), "_", x, ".bam ",
                              "in temp folder."))
                 system(paste0("samtools view -b ", 
                               i,
                               " ", x, " > ",
                               tempfolder, "/", 
                               basename(file_path_sans_ext(i)), "_", x, ".bam"))
                 setwd(tempfolder)
                 sortBam(paste0(basename(file_path_sans_ext(i)), "_", x, ".bam"),
                         paste0(basename(file_path_sans_ext(i)), "_", x))
                 indexBam(paste0(basename(file_path_sans_ext(i)), "_", x, ".bam"))
                 setwd(wd)
               }
               chrbam <- list_files_with_exts(tempfolder, "bam")
               chrbam <- sort(chrbam)
               
               
               for (bam in chrbam) {
                 print(paste0("Identifying variants in ", basename(bam)))
                 bamfl <- BamFile(bam)
                 bpp <- BiocParallel::MulticoreParam(threads)
                 scan <- scanBam(bam)
                 if (length(scan[[1]]$rname) < 1) next
                 chromosome <- as.character(scan[[1]]$rname[[1]])
                 rm(scan)
                 gr <- GRanges(seqinfo(refgen))
                 gr <- gr[seqnames(gr) != chromosome]
                 tally.Param <- TallyVariantsParam(refgen, high_base_quality = 23L,
                                                   mask = gr)
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
                 
                 if (length(snp) > 0)  {
                   vcf <- asVCF(snp)
                   print("Writing to VCF file...")
                   setwd(tempfolder)
                   print(getwd())
                   print(paste0(file_path_sans_ext(basename(bam)), ".vcf"))
                   a <- writeVcf(vcf, paste0(file_path_sans_ext(basename(bam)), ".vcf"), 
                                 index = FALSE)
                   gc()
                 } else {
                   print(paste0("No variants found in ", (basename(bam))))
                 }
                 
                 setwd(wd)
               } 
               
               print("Merging VCF files...")
               vcflist <- list_files_with_exts(tempfolder, "vcf")
               vcfpath <- normalizePath(vcflist)
               match <- VariantAnnotation::readVcf(vcfpath[1])
               print(paste0(gsub(".vcf", "", basename(vcfpath[1])), " OKAY"))
               for (item in vcfpath[-1]) {
                 b <- VariantAnnotation::readVcf(item)
                 match <- BiocGenerics::rbind(match, b)
                 print(paste0(gsub(".vcf", "", basename(item)), " OKAY"))
               }
               rsIDs <- rownames(match)
               rownames(match) <- gsub("chr?(.*)", ".", rsIDs)
               writeVcf(match, paste0(outputdir, "/",
                                      basename(file_path_sans_ext(i)),
                                      "_SNP", ".vcf"), 
                        index = TRUE)
               file.remove(list.files(normalizePath(tempfolder), full.names = TRUE))
               print(paste0("Created ", file_path_sans_ext(basename(i)), "_SNP", ".vcf.bgz"))
               gc()
               setwd(wd)
               append(output, paste0(outputdir, "/", 
                                     file_path_sans_ext(basename(i)),
                                     "_SNP", ".vcf.bgz"))
             }
             file.remove(normalizePath(tempfolder))
           },
           "5"={
             ##Drosophilia melanogaster dm6
             library(BSgenome.Dmelanogaster.UCSC.dm6)
             organism <- BSgenome.Dmelanogaster.UCSC.dm6
             
             ##Selects dm6 as the reference genome
             ##If reference doesn't exist within package directory, create one
             if (dir.exists(paste0(find.package("Exvar"), "/dm6"))) {
               refgen <- GmapGenome(organism,
                                    directory = find.package("Exvar"))
             } else {
               print("Reference genome not found. Creating reference. This might take a while...")
               refgen <- GmapGenome(organism,
                                    directory = find.package("Exvar"),
                                    create = TRUE)
             }
             output <- c()
             for (i in bam) {
               print(paste0("Analysing ", basename(i)))
               if (isFALSE(file.exists(paste0(i, ".bai")))) {
                 print("Indexing bam...")
                 setwd(dirname(i))
                 sortBam(i, file_path_sans_ext(i))
                 indexBam(i)
               }  
               setwd(wd)
               bamfl <- BamFile(i)
               chromosomes <- seqlevels(bamfl)
               print("Splitting BAM into chromosomes...")
               tempfolder <- paste0(dirname(i), "/temp")
               dir.create(tempfolder)
               for (x in chromosomes) {
                 print(paste0("Creating ", 
                              basename(file_path_sans_ext(i)), "_", x, ".bam ",
                              "in temp folder."))
                 system(paste0("samtools view -b ", 
                               i,
                               " ", x, " > ",
                               tempfolder, "/", 
                               basename(file_path_sans_ext(i)), "_", x, ".bam"))
                 setwd(tempfolder)
                 sortBam(paste0(basename(file_path_sans_ext(i)), "_", x, ".bam"),
                         paste0(basename(file_path_sans_ext(i)), "_", x))
                 indexBam(paste0(basename(file_path_sans_ext(i)), "_", x, ".bam"))
                 setwd(wd)
               }
               chrbam <- list_files_with_exts(tempfolder, "bam")
               chrbam <- sort(chrbam)
               
               
               for (bam in chrbam) {
                 print(paste0("Identifying variants in ", basename(bam)))
                 bamfl <- BamFile(bam)
                 bpp <- BiocParallel::MulticoreParam(threads)
                 scan <- scanBam(bam)
                 if (length(scan[[1]]$rname) < 1) next
                 chromosome <- as.character(scan[[1]]$rname[[1]])
                 rm(scan)
                 gr <- GRanges(seqinfo(refgen))
                 gr <- gr[seqnames(gr) != chromosome]
                 tally.Param <- TallyVariantsParam(refgen, high_base_quality = 23L,
                                                   mask = gr)
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
                 
                 if (length(snp) > 0)  {
                   vcf <- asVCF(snp)
                   print("Writing to VCF file...")
                   setwd(tempfolder)
                   print(getwd())
                   print(paste0(file_path_sans_ext(basename(bam)), ".vcf"))
                   a <- writeVcf(vcf, paste0(file_path_sans_ext(basename(bam)), ".vcf"), 
                                 index = FALSE)
                   gc()
                 } else {
                   print(paste0("No variants found in ", (basename(bam))))
                 }
                 
                 setwd(wd)
               } 
               
               print("Merging VCF files...")
               vcflist <- list_files_with_exts(tempfolder, "vcf")
               vcfpath <- normalizePath(vcflist)
               match <- VariantAnnotation::readVcf(vcfpath[1])
               print(paste0(gsub(".vcf", "", basename(vcfpath[1])), " OKAY"))
               for (item in vcfpath[-1]) {
                 b <- VariantAnnotation::readVcf(item)
                 match <- BiocGenerics::rbind(match, b)
                 print(paste0(gsub(".vcf", "", basename(item)), " OKAY"))
               }
               rsIDs <- rownames(match)
               rownames(match) <- gsub("chr?(.*)", ".", rsIDs)
               writeVcf(match, paste0(outputdir, "/",
                                      basename(file_path_sans_ext(i)),
                                      "_SNP", ".vcf"), 
                        index = TRUE)
               file.remove(list.files(normalizePath(tempfolder), full.names = TRUE))
               print(paste0("Created ", file_path_sans_ext(basename(i)), "_SNP", ".vcf.bgz"))
               gc()
               setwd(wd)
               append(output, paste0(outputdir, "/", 
                                     file_path_sans_ext(basename(i)),
                                     "_SNP", ".vcf.bgz"))
             }
             file.remove(normalizePath(tempfolder))
           },
           "6"={
             ##Danio rerio danRer11
             library(BSgenome.Drerio.UCSC.danRer11)
             organism <- BSgenome.Drerio.UCSC.danRer11
             
             ##Selects danRer11 as the reference genome
             ##If reference doesn't exist within package directory, create one
             if (dir.exists(paste0(find.package("Exvar"), "/danRer11"))) {
               refgen <- GmapGenome(organism,
                                    directory = find.package("Exvar"))
             } else {
               print("Reference genome not found. Creating reference. This might take a while...")
               refgen <- GmapGenome(organism,
                                    directory = find.package("Exvar"),
                                    create = TRUE)
             }
             output <- c()
             for (i in bam) {
               print(paste0("Analysing ", basename(i)))
               if (isFALSE(file.exists(paste0(i, ".bai")))) {
                 print("Indexing bam...")
                 setwd(dirname(i))
                 sortBam(i, file_path_sans_ext(i))
                 indexBam(i)
               }  
               setwd(wd)
               bamfl <- BamFile(i)
               chromosomes <- seqlevels(bamfl)
               print("Splitting BAM into chromosomes...")
               tempfolder <- paste0(dirname(i), "/temp")
               dir.create(tempfolder)
               for (x in chromosomes) {
                 print(paste0("Creating ", 
                              basename(file_path_sans_ext(i)), "_", x, ".bam ",
                              "in temp folder."))
                 system(paste0("samtools view -b ", 
                               i,
                               " ", x, " > ",
                               tempfolder, "/", 
                               basename(file_path_sans_ext(i)), "_", x, ".bam"))
                 setwd(tempfolder)
                 sortBam(paste0(basename(file_path_sans_ext(i)), "_", x, ".bam"),
                         paste0(basename(file_path_sans_ext(i)), "_", x))
                 indexBam(paste0(basename(file_path_sans_ext(i)), "_", x, ".bam"))
                 setwd(wd)
               }
               chrbam <- list_files_with_exts(tempfolder, "bam")
               chrbam <- sort(chrbam)
               
               
               for (bam in chrbam) {
                 print(paste0("Identifying variants in ", basename(bam)))
                 bamfl <- BamFile(bam)
                 bpp <- BiocParallel::MulticoreParam(threads)
                 scan <- scanBam(bam)
                 if (length(scan[[1]]$rname) < 1) next
                 chromosome <- as.character(scan[[1]]$rname[[1]])
                 rm(scan)
                 gr <- GRanges(seqinfo(refgen))
                 gr <- gr[seqnames(gr) != chromosome]
                 tally.Param <- TallyVariantsParam(refgen, high_base_quality = 23L,
                                                   mask = gr)
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
                 
                 if (length(snp) > 0)  {
                   vcf <- asVCF(snp)
                   print("Writing to VCF file...")
                   setwd(tempfolder)
                   print(getwd())
                   print(paste0(file_path_sans_ext(basename(bam)), ".vcf"))
                   a <- writeVcf(vcf, paste0(file_path_sans_ext(basename(bam)), ".vcf"), 
                                 index = FALSE)
                   gc()
                 } else {
                   print(paste0("No variants found in ", (basename(bam))))
                 }
                 
                 setwd(wd)
               } 
               
               print("Merging VCF files...")
               vcflist <- list_files_with_exts(tempfolder, "vcf")
               vcfpath <- normalizePath(vcflist)
               match <- VariantAnnotation::readVcf(vcfpath[1])
               print(paste0(gsub(".vcf", "", basename(vcfpath[1])), " OKAY"))
               for (item in vcfpath[-1]) {
                 b <- VariantAnnotation::readVcf(item)
                 match <- BiocGenerics::rbind(match, b)
                 print(paste0(gsub(".vcf", "", basename(item)), " OKAY"))
               }
               rsIDs <- rownames(match)
               rownames(match) <- gsub("chr?(.*)", ".", rsIDs)
               writeVcf(match, paste0(outputdir, "/",
                                      basename(file_path_sans_ext(i)),
                                      "_SNP", ".vcf"), 
                        index = TRUE)
               file.remove(list.files(normalizePath(tempfolder), full.names = TRUE))
               print(paste0("Created ", file_path_sans_ext(basename(i)), "_SNP", ".vcf.bgz"))
               gc()
               setwd(wd)
               append(output, paste0(outputdir, "/", 
                                     file_path_sans_ext(basename(i)),
                                     "_SNP", ".vcf.bgz"))
             }
             file.remove(normalizePath(tempfolder))
           },
           "7"={
             ##Rattus norvegicus rn5
             library(BSgenome.Rnorvegicus.UCSC.rn5)
             organism <- BSgenome.Rnorvegicus.UCSC.rn5
             
             ##Selects danRer11 as the reference genome
             ##If reference doesn't exist within package directory, create one
             if (dir.exists(paste0(find.package("Exvar"), "/rn5"))) {
               refgen <- GmapGenome(organism,
                                    directory = find.package("Exvar"))
             } else {
               print("Reference genome not found. Creating reference. This might take a while...")
               refgen <- GmapGenome(organism,
                                    directory = find.package("Exvar"),
                                    create = TRUE)
             }
             output <- c()
             for (i in bam) {
               print(paste0("Analysing ", basename(i)))
               if (isFALSE(file.exists(paste0(i, ".bai")))) {
                 print("Indexing bam...")
                 setwd(dirname(i))
                 sortBam(i, file_path_sans_ext(i))
                 indexBam(i)
               }  
               setwd(wd)
               bamfl <- BamFile(i)
               chromosomes <- seqlevels(bamfl)
               print("Splitting BAM into chromosomes...")
               tempfolder <- paste0(dirname(i), "/temp")
               dir.create(tempfolder)
               for (x in chromosomes) {
                 print(paste0("Creating ", 
                              basename(file_path_sans_ext(i)), "_", x, ".bam ",
                              "in temp folder."))
                 system(paste0("samtools view -b ", 
                               i,
                               " ", x, " > ",
                               tempfolder, "/", 
                               basename(file_path_sans_ext(i)), "_", x, ".bam"))
                 setwd(tempfolder)
                 sortBam(paste0(basename(file_path_sans_ext(i)), "_", x, ".bam"),
                         paste0(basename(file_path_sans_ext(i)), "_", x))
                 indexBam(paste0(basename(file_path_sans_ext(i)), "_", x, ".bam"))
                 setwd(wd)
               }
               chrbam <- list_files_with_exts(tempfolder, "bam")
               chrbam <- sort(chrbam)
               
               
               for (bam in chrbam) {
                 print(paste0("Identifying variants in ", basename(bam)))
                 bamfl <- BamFile(bam)
                 bpp <- BiocParallel::MulticoreParam(threads)
                 scan <- scanBam(bam)
                 if (length(scan[[1]]$rname) < 1) next
                 chromosome <- as.character(scan[[1]]$rname[[1]])
                 rm(scan)
                 gr <- GRanges(seqinfo(refgen))
                 gr <- gr[seqnames(gr) != chromosome]
                 tally.Param <- TallyVariantsParam(refgen, high_base_quality = 23L,
                                                   mask = gr)
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
                 
                 if (length(snp) > 0)  {
                   vcf <- asVCF(snp)
                   print("Writing to VCF file...")
                   setwd(tempfolder)
                   print(getwd())
                   print(paste0(file_path_sans_ext(basename(bam)), ".vcf"))
                   a <- writeVcf(vcf, paste0(file_path_sans_ext(basename(bam)), ".vcf"), 
                                 index = FALSE)
                   gc()
                 } else {
                   print(paste0("No variants found in ", (basename(bam))))
                 }
                 
                 setwd(wd)
               } 
               
               print("Merging VCF files...")
               vcflist <- list_files_with_exts(tempfolder, "vcf")
               vcfpath <- normalizePath(vcflist)
               match <- VariantAnnotation::readVcf(vcfpath[1])
               print(paste0(gsub(".vcf", "", basename(vcfpath[1])), " OKAY"))
               for (item in vcfpath[-1]) {
                 b <- VariantAnnotation::readVcf(item)
                 match <- BiocGenerics::rbind(match, b)
                 print(paste0(gsub(".vcf", "", basename(item)), " OKAY"))
               }
               rsIDs <- rownames(match)
               rownames(match) <- gsub("chr?(.*)", ".", rsIDs)
               writeVcf(match, paste0(outputdir, "/",
                                      basename(file_path_sans_ext(i)),
                                      "_SNP", ".vcf"), 
                        index = TRUE)
               file.remove(list.files(normalizePath(tempfolder), full.names = TRUE))
               print(paste0("Created ", file_path_sans_ext(basename(i)), "_SNP", ".vcf.bgz"))
               gc()
               setwd(wd)
               append(output, paste0(outputdir, "/", 
                                     file_path_sans_ext(basename(i)),
                                     "_SNP", ".vcf.bgz"))
             }
             file.remove(normalizePath(tempfolder))
           },
           "8"={
             ##Saccharomyces cerevisiae sacCer3
             library(BSgenome.Scerevisiae.UCSC.sacCer3)
             organism <- BSgenome.Scerevisiae.UCSC.sacCer3
             
             ##Selects sacCer3 as the reference genome
             ##If reference doesn't exist within package directory, create one
             if (dir.exists(paste0(find.package("Exvar"), "/sacCer3"))) {
               refgen <- GmapGenome(organism,
                                    directory = find.package("Exvar"))
             } else {
               print("Reference genome not found. Creating reference. This might take a while...")
               refgen <- GmapGenome(organism,
                                    directory = find.package("Exvar"),
                                    create = TRUE)
             }
             
             for (i in bam) {
               print(paste0("Analysing ", basename(i)))
               if (isFALSE(file.exists(paste0(i, ".bai")))) {
                 print("Indexing bam...")
                 setwd(dirname(i))
                 sortBam(i, file_path_sans_ext(i))
                 indexBam(i)
               }  
               setwd(wd)
               bamfl <- BamFile(i)
               chromosomes <- seqlevels(bamfl)
               print("Splitting BAM into chromosomes...")
               tempfolder <- paste0(dirname(i), "/temp")
               dir.create(tempfolder)
               for (x in chromosomes) {
                 print(paste0("Creating ", 
                              basename(file_path_sans_ext(i)), "_", x, ".bam ",
                              "in temp folder."))
                 system(paste0("samtools view -b ", 
                               i,
                               " ", x, " > ",
                               tempfolder, "/", 
                               basename(file_path_sans_ext(i)), "_", x, ".bam"))
                 setwd(tempfolder)
                 sortBam(paste0(basename(file_path_sans_ext(i)), "_", x, ".bam"),
                         paste0(basename(file_path_sans_ext(i)), "_", x))
                 indexBam(paste0(basename(file_path_sans_ext(i)), "_", x, ".bam"))
                 setwd(wd)
               }
               chrbam <- list_files_with_exts(tempfolder, "bam")
               chrbam <- sort(chrbam)
               
               
               for (bam in chrbam) {
                 print(paste0("Identifying variants in ", basename(bam)))
                 bamfl <- BamFile(bam)
                 bpp <- BiocParallel::MulticoreParam(threads)
                 scan <- scanBam(bam)
                 if (length(scan[[1]]$rname) < 1) next
                 chromosome <- as.character(scan[[1]]$rname[[1]])
                 rm(scan)
                 gr <- GRanges(seqinfo(refgen))
                 gr <- gr[seqnames(gr) != chromosome]
                 tally.Param <- TallyVariantsParam(refgen, high_base_quality = 23L,
                                                   mask = gr)
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
                 
                 if (length(snp) > 0)  {
                   vcf <- asVCF(snp)
                   print("Writing to VCF file...")
                   setwd(tempfolder)
                   print(getwd())
                   print(paste0(file_path_sans_ext(basename(bam)), ".vcf"))
                   a <- writeVcf(vcf, paste0(file_path_sans_ext(basename(bam)), ".vcf"), 
                                 index = FALSE)
                   gc()
                 } else {
                   print(paste0("No variants found in ", (basename(bam))))
                 }
                 
                 setwd(wd)
               } 
               
               print("Merging VCF files...")
               vcflist <- list_files_with_exts(tempfolder, "vcf")
               vcfpath <- normalizePath(vcflist)
               match <- VariantAnnotation::readVcf(vcfpath[1])
               print(paste0(gsub(".vcf", "", basename(vcfpath[1])), " OKAY"))
               for (item in vcfpath[-1]) {
                 b <- VariantAnnotation::readVcf(item)
                 match <- BiocGenerics::rbind(match, b)
                 print(paste0(gsub(".vcf", "", basename(item)), " OKAY"))
               }
               rsIDs <- rownames(match)
               rownames(match) <- gsub("chr?(.*)", ".", rsIDs)
               writeVcf(match, paste0(outputdir, "/",
                                      basename(file_path_sans_ext(i)),
                                      "_SNP", ".vcf"), 
                        index = TRUE)
               file.remove(list.files(normalizePath(tempfolder), full.names = TRUE))
               print(paste0("Created ", file_path_sans_ext(basename(i)), "_SNP", ".vcf.bgz"))
               gc()
               setwd(wd)
               append(output, paste0(outputdir, "/", 
                                     file_path_sans_ext(basename(i)),
                                     "_SNP", ".vcf.bgz"))
             }
             file.remove(normalizePath(tempfolder))
           },
           "9"={
             ##Caenorhabditis elagans
             library(BSgenome.Celegans.UCSC.ce11)
             organism <- BSgenome.Celegans.UCSC.ce11
             
             ##Selects ce11 as the reference genome
             ##If reference doesn't exist within package directory, create one
             if (dir.exists(paste0(find.package("Exvar"), "/ce11"))) {
               refgen <- GmapGenome(organism,
                                    directory = find.package("Exvar"))
             } else {
               print("Reference genome not found. Creating reference. This might take a while...")
               refgen <- GmapGenome(organism,
                                    directory = find.package("Exvar"),
                                    create = TRUE)
             }
             output <- c()
             for (i in bam) {
               print(paste0("Analysing ", basename(i)))
               if (isFALSE(file.exists(paste0(i, ".bai")))) {
                 print("Indexing bam...")
                 setwd(dirname(i))
                 sortBam(i, file_path_sans_ext(i))
                 indexBam(i)
               }  
               setwd(wd)
               bamfl <- BamFile(i)
               chromosomes <- seqlevels(bamfl)
               print("Splitting BAM into chromosomes...")
               tempfolder <- paste0(dirname(i), "/temp")
               dir.create(tempfolder)
               for (x in chromosomes) {
                 print(paste0("Creating ", 
                              basename(file_path_sans_ext(i)), "_", x, ".bam ",
                              "in temp folder."))
                 system(paste0("samtools view -b ", 
                               i,
                               " ", x, " > ",
                               tempfolder, "/", 
                               basename(file_path_sans_ext(i)), "_", x, ".bam"))
                 setwd(tempfolder)
                 sortBam(paste0(basename(file_path_sans_ext(i)), "_", x, ".bam"),
                         paste0(basename(file_path_sans_ext(i)), "_", x))
                 indexBam(paste0(basename(file_path_sans_ext(i)), "_", x, ".bam"))
                 setwd(wd)
               }
               chrbam <- list_files_with_exts(tempfolder, "bam")
               chrbam <- sort(chrbam)
               
               
               for (bam in chrbam) {
                 print(paste0("Identifying variants in ", basename(bam)))
                 bamfl <- BamFile(bam)
                 bpp <- BiocParallel::MulticoreParam(threads)
                 scan <- scanBam(bam)
                 if (length(scan[[1]]$rname) < 1) next
                 chromosome <- as.character(scan[[1]]$rname[[1]])
                 rm(scan)
                 gr <- GRanges(seqinfo(refgen))
                 gr <- gr[seqnames(gr) != chromosome]
                 tally.Param <- TallyVariantsParam(refgen, high_base_quality = 23L,
                                                   mask = gr)
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
                 
                 if (length(snp) > 0)  {
                   vcf <- asVCF(snp)
                   print("Writing to VCF file...")
                   setwd(tempfolder)
                   print(getwd())
                   print(paste0(file_path_sans_ext(basename(bam)), ".vcf"))
                   a <- writeVcf(vcf, paste0(file_path_sans_ext(basename(bam)), ".vcf"), 
                                 index = FALSE)
                   gc()
                 } else {
                   print(paste0("No variants found in ", (basename(bam))))
                 }
                 
                 setwd(wd)
               } 
               
               print("Merging VCF files...")
               vcflist <- list_files_with_exts(tempfolder, "vcf")
               vcfpath <- normalizePath(vcflist)
               match <- VariantAnnotation::readVcf(vcfpath[1])
               print(paste0(gsub(".vcf", "", basename(vcfpath[1])), " OKAY"))
               for (item in vcfpath[-1]) {
                 b <- VariantAnnotation::readVcf(item)
                 match <- BiocGenerics::rbind(match, b)
                 print(paste0(gsub(".vcf", "", basename(item)), " OKAY"))
               }
               rsIDs <- rownames(match)
               rownames(match) <- gsub("chr?(.*)", ".", rsIDs)
               writeVcf(match, paste0(outputdir, "/",
                                      basename(file_path_sans_ext(i)),
                                      "_SNP", ".vcf"), 
                        index = TRUE)
               file.remove(list.files(normalizePath(tempfolder), full.names = TRUE))
               print(paste0("Created ", file_path_sans_ext(basename(i)), "_SNP", ".vcf.bgz"))
               gc()
               setwd(wd)
               append(output, paste0(outputdir, "/", 
                                     file_path_sans_ext(basename(i)),
                                     "_SNP", ".vcf.bgz"))
             }
             file.remove(normalizePath(tempfolder))
           }
    )
  } else if(mode == "n") {
    ##Sets the reference genome that corresponds to the species chosen by the user
    switch(species,
           "1"={
             ##Homo sapiens hg19
             library(BSgenome.Hsapiens.UCSC.hg19)
             library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
             organism <- BSgenome.Hsapiens.UCSC.hg19
             
             ##Selects hg19 as the reference genome
             ##If reference doesn't exist within package directory, create one
             if (dir.exists(paste0(find.package("Exvar"), "/hg19"))) {
               refgen <- GmapGenome(organism,
                                    directory = find.package("Exvar"))
             } else {
               print("Reference genome not found. Creating reference. This might take a while...")
               refgen <- GmapGenome(organism,
                                    directory = find.package("Exvar"),
                                    create = TRUE)
             }
             
             output <- c()
             for(i in bam) {
               print(paste0("Analysing ", basename(i)))
               if (isFALSE(file.exists(paste0(i, ".bai")))) {
                 print("Indexing bam...")
                 setwd(dirname(i))
                 sortBam(i, file_path_sans_ext(i))
                 indexBam(i)
               } 
               tempfolder <- paste0(dirname(i), "/temp")
               dir.create(tempfolder)
               vcflist <- c()
               bpp <- BiocParallel::MulticoreParam(threads)
               chromosomes <- seqlevels(refgen)
               bamfl <- BamFile(i)
               tally.Param <- TallyVariantsParam(refgen, high_base_quality = 23L)
               print("Tallying BAM file...")
               tallies <- tallyVariants(bamfl, tally.Param, BPPARAM = bpp)
               gc()
               calling.filters <- VariantCallingFilters()
               post.filters <- VariantPostFilters()
               print("Calling variants...")
               snp <- callVariants(tallies, calling.filters, post.filters)
               sampleNames(snp) <- file_path_sans_ext(dirname(i))
               mcols(snp) <- NULL
               print("Formatting variant information as VCF...")
               vcf <- asVCF(snp)
               vcflist <- c()
               
               writeVcf(vcf, paste0(tempfolder, "/",
                                    basename(file_path_sans_ext(i)),
                                    "_SNP", ".vcf"),
                        index = FALSE)
               vcf <- readVcf(paste0(tempfolder, "/",
                                     basename(file_path_sans_ext(i)),
                                     "_SNP", ".vcf"))
               
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
               writeVcf(vcf, paste0(outputdir, "/",
                                    basename(file_path_sans_ext(i)),
                                    "_SNP", ".vcf"), 
                        index = TRUE)
               file.remove(list.files(normalizePath(tempfolder), full.names = TRUE))
               print(paste0("Created ", file_path_sans_ext(basename(i)), "_SNP", ".vcf.bgz"))
               file.remove(normalizePath(tempfolder))
               gc()
               setwd(wd)
               append(output, paste0(outputdir, "/", 
                                     file_path_sans_ext(basename(i)),
                                     "_SNP", ".vcf.bgz"))
               
             }
             
           },
           "2"={
             ##Homo sapiens hg38
             library(BSgenome.Hsapiens.UCSC.hg38)
             library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
             organism <- BSgenome.Hsapiens.UCSC.hg38
             
             ##Selects hg38 as the reference genome
             ##If reference doesn't exist within package directory, create one
             if (dir.exists(paste0(find.package("Exvar"), "/hg38"))) {
               refgen <- GmapGenome(organism,
                                    directory = find.package("Exvar"))
             } else {
               print("Reference genome not found. Creating reference. This might take a while...")
               refgen <- GmapGenome(organism,
                                    directory = find.package("Exvar"),
                                    create = TRUE)
             }
             output <- c()
             for(i in bam) {
               print(paste0("Analysing ", basename(i)))
               if (isFALSE(file.exists(paste0(i, ".bai")))) {
                 print("Indexing bam...")
                 setwd(dirname(i))
                 sortBam(i, file_path_sans_ext(i))
                 indexBam(i)
               } 
               tempfolder <- paste0(dirname(i), "/temp")
               dir.create(tempfolder)
               vcflist <- c()
               bpp <- BiocParallel::MulticoreParam(threads)
               chromosomes <- seqlevels(refgen)
               bamfl <- BamFile(i)
               tally.Param <- TallyVariantsParam(refgen, high_base_quality = 23L)
               print("Tallying BAM file...")
               tallies <- tallyVariants(bamfl, tally.Param, BPPARAM = bpp)
               gc()
               calling.filters <- VariantCallingFilters()
               post.filters <- VariantPostFilters()
               print("Calling variants...")
               snp <- callVariants(tallies, calling.filters, post.filters)
               sampleNames(snp) <- file_path_sans_ext(dirname(i))
               mcols(snp) <- NULL
               print("Formatting variant information as VCF...")
               vcf <- asVCF(snp)
               writeVcf(vcf, paste0(tempfolder, "/",
                                    basename(file_path_sans_ext(i)),
                                    "_SNP", ".vcf"))
               
               
               vcflist <- c()
               vcf <- readVcf(paste0(tempfolder, "/",
                                     basename(file_path_sans_ext(i)),
                                     "_SNP", ".vcf"))
               
               ##Due to constraints in memory, rsIDs are obtained on a chromosome by chromosome
               ##basis
               ##The resulting data frame is sorted as per the vcf and the rsIDs are injected
               ##into the vcf before writing it as a file again
               print("Loading dbSNP information...")
               all_snps <- SNPlocs.Hsapiens.dbSNP151.GRCh38
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
               writeVcf(vcf, paste0(outputdir, "/",
                                    basename(file_path_sans_ext(i)),
                                    "_SNP", ".vcf"), 
                        index = TRUE)
               file.remove(list.files(normalizePath(tempfolder), full.names = TRUE))
               print(paste0("Created ", file_path_sans_ext(basename(i)), "_SNP", ".vcf.bgz"))
               file.remove(normalizePath(tempfolder))
               gc()
               setwd(wd)
               append(output, paste0(outputdir, "/", 
                                     file_path_sans_ext(basename(i)),
                                     "_SNP", ".vcf.bgz"))
               
             }
           },
           "3"={
             ##Mus musculus mm10
             library(BSgenome.Mmusculus.UCSC.mm10)
             organism <- BSgenome.Mmusculus.UCSC.mm10
             
             ##Selects mm10 as the reference genome
             ##If reference doesn't exist within package directory, create one
             if (dir.exists(paste0(find.package("Exvar"), "/mm10"))) {
               refgen <- GmapGenome(organism,
                                    directory = find.package("Exvar"))
             } else {
               print("Reference genome not found. Creating reference. This might take a while...")
               refgen <- GmapGenome(organism,
                                    directory = find.package("Exvar"),
                                    create = TRUE)
             }
             output <- c()
             output <- c()
             for(i in bam) {
               print(paste0("Analysing ", basename(i)))
               if (isFALSE(file.exists(paste0(i, ".bai")))) {
                 print("Indexing bam...")
                 setwd(dirname(i))
                 sortBam(i, file_path_sans_ext(i))
                 indexBam(i)
               } 
               tempfolder <- paste0(dirname(i), "/temp")
               dir.create(tempfolder)
               vcflist <- c()
               bpp <- BiocParallel::MulticoreParam(threads)
               chromosomes <- seqlevels(refgen)
               bamfl <- BamFile(i)
               tally.Param <- TallyVariantsParam(refgen, high_base_quality = 23L)
               print("Tallying BAM file...")
               tallies <- tallyVariants(bamfl, tally.Param, BPPARAM = bpp)
               gc()
               calling.filters <- VariantCallingFilters()
               post.filters <- VariantPostFilters()
               print("Calling variants...")
               snp <- callVariants(tallies, calling.filters, post.filters)
               sampleNames(snp) <- file_path_sans_ext(dirname(i))
               mcols(snp) <- NULL
               print("Formatting variant information as VCF...")
               vcf <- asVCF(snp)
               print("Writing to VCF file...")
               setwd(outputdir)
               writeVcf(vcf, paste0(outputdir, "/",
                                    basename(file_path_sans_ext(i)),
                                    "_SNP", ".vcf"), 
                        index = TRUE)
               file.remove(list.files(normalizePath(tempfolder), full.names = TRUE))
               print(paste0("Created ", file_path_sans_ext(basename(i)), "_SNP", ".vcf.bgz"))
               file.remove(normalizePath(tempfolder))
               gc()
               setwd(wd)
               append(output, paste0(outputdir, "/", 
                                     file_path_sans_ext(basename(i)),
                                     "_SNP", ".vcf.bgz"))
             }
           },
           "4"={
             ##Arabidopsis thaliana TAIR9
             library(BSgenome.Athaliana.TAIR.TAIR9)
             organism <- BSgenome.Athaliana.TAIR.TAIR9
             
             ##Selects hg19 as the reference genome
             ##If reference doesn't exist within package directory, create one
             if (dir.exists(paste0(find.package("Exvar"), "/TAIR9"))) {
               refgen <- GmapGenome(organism,
                                    directory = find.package("Exvar"))
             } else {
               print("Reference genome not found. Creating reference. This might take a while...")
               refgen <- GmapGenome(organism,
                                    directory = find.package("Exvar"),
                                    create = TRUE)
             }
             output <- c()
             output <- c()
             for(i in bam) {
               print(paste0("Analysing ", basename(i)))
               if (isFALSE(file.exists(paste0(i, ".bai")))) {
                 print("Indexing bam...")
                 setwd(dirname(i))
                 sortBam(i, file_path_sans_ext(i))
                 indexBam(i)
               } 
               tempfolder <- paste0(dirname(i), "/temp")
               dir.create(tempfolder)
               vcflist <- c()
               bpp <- BiocParallel::MulticoreParam(threads)
               chromosomes <- seqlevels(refgen)
               bamfl <- BamFile(i)
               tally.Param <- TallyVariantsParam(refgen, high_base_quality = 23L)
               print("Tallying BAM file...")
               tallies <- tallyVariants(bamfl, tally.Param, BPPARAM = bpp)
               gc()
               calling.filters <- VariantCallingFilters()
               post.filters <- VariantPostFilters()
               print("Calling variants...")
               snp <- callVariants(tallies, calling.filters, post.filters)
               sampleNames(snp) <- file_path_sans_ext(dirname(i))
               mcols(snp) <- NULL
               print("Formatting variant information as VCF...")
               vcf <- asVCF(snp)
               print("Writing to VCF file...")
               setwd(outputdir)
               writeVcf(vcf, paste0(outputdir, "/",
                                    basename(file_path_sans_ext(i)),
                                    "_SNP", ".vcf"), 
                        index = TRUE)
               file.remove(list.files(normalizePath(tempfolder), full.names = TRUE))
               print(paste0("Created ", file_path_sans_ext(basename(i)), "_SNP", ".vcf.bgz"))
               file.remove(normalizePath(tempfolder))
               gc()
               setwd(wd)
               append(output, paste0(outputdir, "/", 
                                     file_path_sans_ext(basename(i)),
                                     "_SNP", ".vcf.bgz"))
             }
           },
           "5"={
             ##Drosophilia melanogaster dm6
             library(BSgenome.Dmelanogaster.UCSC.dm6)
             organism <- BSgenome.Dmelanogaster.UCSC.dm6
             
             ##Selects dm6 as the reference genome
             ##If reference doesn't exist within package directory, create one
             if (dir.exists(paste0(find.package("Exvar"), "/dm6"))) {
               refgen <- GmapGenome(organism,
                                    directory = find.package("Exvar"))
             } else {
               print("Reference genome not found. Creating reference. This might take a while...")
               refgen <- GmapGenome(organism,
                                    directory = find.package("Exvar"),
                                    create = TRUE)
             }
             output <- c()
             output <- c()
             for(i in bam) {
               print(paste0("Analysing ", basename(i)))
               if (isFALSE(file.exists(paste0(i, ".bai")))) {
                 print("Indexing bam...")
                 setwd(dirname(i))
                 sortBam(i, file_path_sans_ext(i))
                 indexBam(i)
               } 
               tempfolder <- paste0(dirname(i), "/temp")
               dir.create(tempfolder)
               vcflist <- c()
               bpp <- BiocParallel::MulticoreParam(threads)
               chromosomes <- seqlevels(refgen)
               bamfl <- BamFile(i)
               tally.Param <- TallyVariantsParam(refgen, high_base_quality = 23L)
               print("Tallying BAM file...")
               tallies <- tallyVariants(bamfl, tally.Param, BPPARAM = bpp)
               gc()
               calling.filters <- VariantCallingFilters()
               post.filters <- VariantPostFilters()
               print("Calling variants...")
               snp <- callVariants(tallies, calling.filters, post.filters)
               sampleNames(snp) <- file_path_sans_ext(dirname(i))
               mcols(snp) <- NULL
               print("Formatting variant information as VCF...")
               vcf <- asVCF(snp)
               print("Writing to VCF file...")
               setwd(outputdir)
               writeVcf(vcf, paste0(outputdir, "/",
                                    basename(file_path_sans_ext(i)),
                                    "_SNP", ".vcf"), 
                        index = TRUE)
               file.remove(list.files(normalizePath(tempfolder), full.names = TRUE))
               print(paste0("Created ", file_path_sans_ext(basename(i)), "_SNP", ".vcf.bgz"))
               file.remove(normalizePath(tempfolder))
               gc()
               setwd(wd)
               append(output, paste0(outputdir, "/", 
                                     file_path_sans_ext(basename(i)),
                                     "_SNP", ".vcf.bgz"))
             }
           },
           "6"={
             ##Danio rerio danRer11
             library(BSgenome.Drerio.UCSC.danRer11)
             organism <- BSgenome.Drerio.UCSC.danRer11
             
             ##Selects danRer11 as the reference genome
             ##If reference doesn't exist within package directory, create one
             if (dir.exists(paste0(find.package("Exvar"), "/danRer11"))) {
               refgen <- GmapGenome(organism,
                                    directory = find.package("Exvar"))
             } else {
               print("Reference genome not found. Creating reference. This might take a while...")
               refgen <- GmapGenome(organism,
                                    directory = find.package("Exvar"),
                                    create = TRUE)
             }
             output <- c()
             for(i in bam) {
               print(paste0("Analysing ", basename(i)))
               if (isFALSE(file.exists(paste0(i, ".bai")))) {
                 print("Indexing bam...")
                 setwd(dirname(i))
                 sortBam(i, file_path_sans_ext(i))
                 indexBam(i)
               } 
               tempfolder <- paste0(dirname(i), "/temp")
               dir.create(tempfolder)
               vcflist <- c()
               bpp <- BiocParallel::MulticoreParam(threads)
               chromosomes <- seqlevels(refgen)
               bamfl <- BamFile(i)
               tally.Param <- TallyVariantsParam(refgen, high_base_quality = 23L)
               print("Tallying BAM file...")
               tallies <- tallyVariants(bamfl, tally.Param, BPPARAM = bpp)
               gc()
               calling.filters <- VariantCallingFilters()
               post.filters <- VariantPostFilters()
               print("Calling variants...")
               snp <- callVariants(tallies, calling.filters, post.filters)
               sampleNames(snp) <- file_path_sans_ext(dirname(i))
               mcols(snp) <- NULL
               print("Formatting variant information as VCF...")
               vcf <- asVCF(snp)
               print("Writing to VCF file...")
               setwd(outputdir)
               writeVcf(vcf, paste0(outputdir, "/",
                                    basename(file_path_sans_ext(i)),
                                    "_SNP", ".vcf"), 
                        index = TRUE)
               file.remove(list.files(normalizePath(tempfolder), full.names = TRUE))
               print(paste0("Created ", file_path_sans_ext(basename(i)), "_SNP", ".vcf.bgz"))
               file.remove(normalizePath(tempfolder))
               gc()
               setwd(wd)
               append(output, paste0(outputdir, "/", 
                                     file_path_sans_ext(basename(i)),
                                     "_SNP", ".vcf.bgz"))
             }
           },
           "7"={
             ##Rattus norvegicus rn5
             library(BSgenome.Rnorvegicus.UCSC.rn5)
             organism <- BSgenome.Rnorvegicus.UCSC.rn5
             
             ##Selects danRer11 as the reference genome
             ##If reference doesn't exist within package directory, create one
             if (dir.exists(paste0(find.package("Exvar"), "/rn5"))) {
               refgen <- GmapGenome(organism,
                                    directory = find.package("Exvar"))
             } else {
               print("Reference genome not found. Creating reference. This might take a while...")
               refgen <- GmapGenome(organism,
                                    directory = find.package("Exvar"),
                                    create = TRUE)
             }
             output <- c()
             for(i in bam) {
               print(paste0("Analysing ", basename(i)))
               if (isFALSE(file.exists(paste0(i, ".bai")))) {
                 print("Indexing bam...")
                 setwd(dirname(i))
                 sortBam(i, file_path_sans_ext(i))
                 indexBam(i)
               } 
               tempfolder <- paste0(dirname(i), "/temp")
               dir.create(tempfolder)
               vcflist <- c()
               bpp <- BiocParallel::MulticoreParam(threads)
               chromosomes <- seqlevels(refgen)
               bamfl <- BamFile(i)
               tally.Param <- TallyVariantsParam(refgen, high_base_quality = 23L)
               print("Tallying BAM file...")
               tallies <- tallyVariants(bamfl, tally.Param, BPPARAM = bpp)
               gc()
               calling.filters <- VariantCallingFilters()
               post.filters <- VariantPostFilters()
               print("Calling variants...")
               snp <- callVariants(tallies, calling.filters, post.filters)
               sampleNames(snp) <- file_path_sans_ext(dirname(i))
               mcols(snp) <- NULL
               print("Formatting variant information as VCF...")
               vcf <- asVCF(snp)
               print("Writing to VCF file...")
               setwd(outputdir)
               writeVcf(vcf, paste0(outputdir, "/",
                                    basename(file_path_sans_ext(i)),
                                    "_SNP", ".vcf"), 
                        index = TRUE)
               file.remove(list.files(normalizePath(tempfolder), full.names = TRUE))
               print(paste0("Created ", file_path_sans_ext(basename(i)), "_SNP", ".vcf.bgz"))
               file.remove(normalizePath(tempfolder))
               gc()
               setwd(wd)
               append(output, paste0(outputdir, "/", 
                                     file_path_sans_ext(basename(i)),
                                     "_SNP", ".vcf.bgz"))
             }
           },
           "8"={
             ##Saccharomyces cerevisiae sacCer3
             library(BSgenome.Scerevisiae.UCSC.sacCer3)
             organism <- BSgenome.Scerevisiae.UCSC.sacCer3
             
             ##Selects sacCer3 as the reference genome
             ##If reference doesn't exist within package directory, create one
             if (dir.exists(paste0(find.package("Exvar"), "/sacCer3"))) {
               refgen <- GmapGenome(organism,
                                    directory = find.package("Exvar"))
             } else {
               print("Reference genome not found. Creating reference. This might take a while...")
               refgen <- GmapGenome(organism,
                                    directory = find.package("Exvar"),
                                    create = TRUE)
             }
             ouput <- c()
             output <- c()
             for(i in bam) {
               print(paste0("Analysing ", basename(i)))
               if (isFALSE(file.exists(paste0(i, ".bai")))) {
                 print("Indexing bam...")
                 setwd(dirname(i))
                 sortBam(i, file_path_sans_ext(i))
                 indexBam(i)
               } 
               tempfolder <- paste0(dirname(i), "/temp")
               dir.create(tempfolder)
               vcflist <- c()
               bpp <- BiocParallel::MulticoreParam(threads)
               chromosomes <- seqlevels(refgen)
               bamfl <- BamFile(i)
               tally.Param <- TallyVariantsParam(refgen, high_base_quality = 23L)
               print("Tallying BAM file...")
               tallies <- tallyVariants(bamfl, tally.Param, BPPARAM = bpp)
               gc()
               calling.filters <- VariantCallingFilters()
               post.filters <- VariantPostFilters()
               print("Calling variants...")
               snp <- callVariants(tallies, calling.filters, post.filters)
               sampleNames(snp) <- file_path_sans_ext(dirname(i))
               mcols(snp) <- NULL
               print("Formatting variant information as VCF...")
               vcf <- asVCF(snp)
               print("Writing to VCF file...")
               setwd(outputdir)
               writeVcf(vcf, paste0(outputdir, "/",
                                    basename(file_path_sans_ext(i)),
                                    "_SNP", ".vcf"), 
                        index = TRUE)
               file.remove(list.files(normalizePath(tempfolder), full.names = TRUE))
               print(paste0("Created ", file_path_sans_ext(basename(i)), "_SNP", ".vcf.bgz"))
               file.remove(normalizePath(tempfolder))
               gc()
               setwd(wd)
               append(output, paste0(outputdir, "/", 
                                     file_path_sans_ext(basename(i)),
                                     "_SNP", ".vcf.bgz"))
             }
           },
           "9"={
             ##Caenorhabditis elagans
             library(BSgenome.Celegans.UCSC.ce11)
             organism <- BSgenome.Celegans.UCSC.ce11
             
             ##Selects ce11 as the reference genome
             ##If reference doesn't exist within package directory, create one
             if (dir.exists(paste0(find.package("Exvar"), "/ce11"))) {
               refgen <- GmapGenome(organism,
                                    directory = find.package("Exvar"))
             } else {
               print("Reference genome not found. Creating reference. This might take a while...")
               refgen <- GmapGenome(organism,
                                    directory = find.package("Exvar"),
                                    create = TRUE)
             }
             output <- c()
             for(i in bam) {
               print(paste0("Analysing ", basename(i)))
               if (isFALSE(file.exists(paste0(i, ".bai")))) {
                 print("Indexing bam...")
                 setwd(dirname(i))
                 sortBam(i, file_path_sans_ext(i))
                 indexBam(i)
               } 
               tempfolder <- paste0(dirname(i), "/temp")
               dir.create(tempfolder)
               vcflist <- c()
               bpp <- BiocParallel::MulticoreParam(threads)
               chromosomes <- seqlevels(refgen)
               bamfl <- BamFile(i)
               tally.Param <- TallyVariantsParam(refgen, high_base_quality = 23L)
               print("Tallying BAM file...")
               tallies <- tallyVariants(bamfl, tally.Param, BPPARAM = bpp)
               gc()
               calling.filters <- VariantCallingFilters()
               post.filters <- VariantPostFilters()
               print("Calling variants...")
               snp <- callVariants(tallies, calling.filters, post.filters)
               sampleNames(snp) <- file_path_sans_ext(dirname(i))
               mcols(snp) <- NULL
               print("Formatting variant information as VCF...")
               vcf <- asVCF(snp)
               print("Writing to VCF file...")
               setwd(outputdir)
               writeVcf(vcf, paste0(outputdir, "/",
                                    basename(file_path_sans_ext(i)),
                                    "_SNP", ".vcf"), 
                        index = TRUE)
               file.remove(list.files(normalizePath(tempfolder), full.names = TRUE))
               print(paste0("Created ", file_path_sans_ext(basename(i)), "_SNP", ".vcf.bgz"))
               file.remove(normalizePath(tempfolder))
               gc()
               setwd(wd)
               append(output, paste0(outputdir, "/", 
                                     file_path_sans_ext(basename(i)),
                                     "_SNP", ".vcf.bgz"))
             }
           }
    )
  }
  
  
  setwd(wd)
  return(output)
}



