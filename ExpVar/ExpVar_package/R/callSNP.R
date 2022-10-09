#' Call single nucleotide polymorphism variants
#' 
#' This function calls SNP variants from BAM files. The results are formatted 
#' into a VCF file and the ID column is populated with dbSNP IDs.
#' 
#' @param bam A list of paths to BAM files
#' @param genome A BSgenome object, a GmapGenome object, or a character string indicating the reference genome eg. "hg19"
#' @param genomedir The directory containing the reference genome or, if genome is a character string, the parent directory of the reference genome directory.
#' @param SNPlocs An SNPlocs object containing dbSNP IDs.
#' @param threads The number of cores to use.
#' @param outputdir The output directory for the VCF file.
#' @return A list of file paths to the VCF files.
#' @export
callSNP <- function(bam,
                    genome,
                    genomedir,
                    SNPlocs,
                    threads = 4L,
                    outputdir = getwd()) {
  
  wd <- getwd()
  if (class(genome) != "GmapGenome") {
    if(isFALSE(exists("genome"))) {
      if (exists("genomedir")) {
        if (dir.exists(genomedir)) {
          refgen <- GmapGenome(BSgenome.Hsapiens.UCSC.hg19, 
                               directory = dirname(genomedir),
                               name = basename(genomedir))
        }
      }
    } else {
      if (isFALSE(exists("genomedir"))) {
        genomedir <- getwd()
      }
      if (class(genome) == "character") {
        if (dir.exists(genome)) {
          refgen <- GmapGenome(BSgenome.Hsapiens.UCSC.hg19, genomedir, 
                               name = genome)
        } 
      } else {
        if (dir.exists(paste0(genomdir, "/", metadata(genome)$genome))) {
          refgen <- GmapGenome(genome, genomedir, create = FALSE) 
        } else {
          outcome <- readline("Reference genome doesn't exist in directory. Create one? [y/n]")
          if (outcome == "y") {
            refgen <- GmapGenome(genome, genomedir, create = TRUE)
          } else {
            stop("Cannot proceed without reference genome.")
          }
        }
      }
    }
  } else {
    refgen <- genome
  }
  vcflist <- c()
  bamfl <- BamFile(bam)
  bpp <- BiocParallel::MulticoreParam(threads)
  tally.Param <- TallyVariantsParam(refgen, high_base_quality = 23L)
  for (i in bamfl) {
    tallies <- tallyVariants(i, tally.Param, BPPARAM = bpp)
    calling.filters <- VariantCallingFilters()
    post.filters <- VariantPostFilters()
    snp <- callVariants(tallies, calling.filters, post.filters)
    sampleNames(snp) <- file_path_sans_ext(basename(path(i)))
    mcols(snp) <- NULL
    print("Formatting variant information as VCF...")
    vcf <- asVCF(snp)
    writeVcf(vcf, paste0(file_path_sans_ext(basename(path(i))), ".vcf"), 
             index = TRUE)
    vcf <- readVcf(paste0(file_path_sans_ext(basename(path(i))), ".vcf.bgz"))
    
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
    a <- writeVcf(vcf, paste0(file_path_sans_ext(basename(path(i))), ".vcf"), 
             index = TRUE)
    a <- paste0(outputdir, "/", a)
    vcflist <- append(vcflist, a)
  }
  setwd(wd)
  return(vcflist)
}





