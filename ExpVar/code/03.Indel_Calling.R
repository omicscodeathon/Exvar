library("gmapR")
library("BSgenome.Hsapiens.UCSC.hg19")
library("Rsamtools")
library("VariantTools")
library("VariantAnnotation")
library("BiocParallel")
library("XtraSNPlocs.Hsapiens.dbSNP144.GRCh37")
library("dplyr")

##Variant calling requires a reference genome so one is created
chr <- 1:22
chrxym <- c("X", "Y", "MT")
chr <- append(chr, chrxym)
mainChromosomes <- paste0("chr", chr)

refSeq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, mainChromosomes)
refgen <- GmapGenome(refSeq, create = TRUE, name = "hg19",
                     directory = getwd())
rm(refSeq)

##Variants are called from a bam file and outputted to a vcf
##The output vcf is a read into memory again so that the rsIDs from dbSNP can
##be integrated into the ID column
##Number of CPUs needs to be a factor of the number of chromosomes
bamfl <- BamFile("Insert bam path here")
print("Tallying variants...")
bpp <- BiocParallel::MulticoreParam(length(seqnames(refgen)))
tally.Param <- TallyVariantsParam(refgen, high_base_quality = 23L, 
                                  indels = TRUE)
tallies <- tallyVariants(bamfl, tally.Param, BPPARAM = bpp)
calling.filters <- VariantCallingFilters()
post.filters <- VariantPostFilters()
indels <- callVariants(tallies, calling.filters, post.filters)
sampleNames(indels) <- file_path_sans_ext(basename(bamfl))
mcols(indels) <- NULL
print("Formatting variant information as VCF...")
vcf <- asVCF(indels)
writeVcf(vcf, 
         paste0(file_path_sans_ext(basename(bamfl)), "_indels.vcf"), 
         index = TRUE)
vcf <- readVcf(paste0(file_path_sans_ext(basename(bamfl)), "_indels.vcf.bgz"))

##SNPs are discarded
##Due to constraints in memory, rsIDs are obtained on a chromosome by chromosome
##basis
##The resulting data frame is sorted as per the vcf and the rsIDs are injected
##into the vcf before writing it as a file again
print("Loading dbSNP information...")
all_snps <- XtraSNPlocs.Hsapiens.dbSNP144.GRCh37
vcfID <- data.frame()
vcf <- vcf[isIndel(vcf)]

print("Finding rsIDs...")
for (x in mainChromosomes) {
  vcf_chrom <- vcf[grepl(names(vcf),
                         pattern = paste0(x, ":"))]
  rd_chr <- rowRanges(vcf_chrom)
  tar_chr <- as.vector(seqnames(rd_chr)@values)
  tar_chr <- gsub("chr", "ch", tar_chr)
  my_snps <- snpsBySeqname(all_snps, c(tar_chr))
  snp_ID <- data.frame(posIDX = paste0(x, 
                                       ":", 
                                       start(my_snps),
                                       "-",
                                       end(my_snps)), 
                       rsID = my_snps$RefSNP_id, 
                       stringsAsFactors = FALSE)
  matV1 <- data.frame(Variant = names(rd_chr), stringsAsFactors = FALSE)
  matV1$chromosome <- gsub("(.*):(.*)_(.*)/(.*)", "\\1", matV1$Variant)
  matV1$ref_allele <- gsub("(.*):(.*)_(.*)/(.*)", "\\3", matV1$Variant)
  matV1$alt_allele <- gsub("(.*):(.*)_(.*)/(.*)", "\\4", matV1$Variant)
  matV1$start <- gsub("(.*):(.*)_(.*)/(.*)", "\\2", matV1$Variant)
  matV1$start <- as.integer(matV1$start)
  length <- 1:length(matV1$start)
  for (i in length) {
    matV1$end[i] <- 
      matV1$start[i] + (nchar(matV1$ref_allele[i]) - nchar(matV1$alt_allele[i]))
  }
  matV1$end[matV1$end < matV1$start] <- matV1$start[matV1$end < matV1$start] - 1
  matV1$posIDX <- gsub("(.*)_(.*)", "\\1", paste0(matV1$chromosome, 
                                                  ":", 
                                                  matV1$start,
                                                  "-",
                                                  matV1$end))
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
vcfID$chr[vcfID$chromosome == "chrMT"] <- 25L
dd <- vcfID[order(vcfID$chr, as.integer(vcfID$start)),]
dd$rsID[is.na(dd$rsID)] <- "."
names(vcf) <- dd$rsID

print("Writing to VCF file...")
writeVcf(vcf, paste0(file_path_sans_ext(basename(bamfl)), "_indels.vcf"), 
         index = TRUE)