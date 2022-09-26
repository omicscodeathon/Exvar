library("gmapR")
library("BSgenome.Hsapiens.UCSC.hg19")
library("Rsamtools")
library("VariantTools")
library("VariantAnnotation")
library("BiocParallel")
library("SNPlocs.Hsapiens.dbSNP144.GRCh37")

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
bamfl <- BamFile("Insert bam path here")
bpp <- BiocParallel::MulticoreParam(28L)
tally.Param <- TallyVariantsParam(refgen, high_base_quality = 23L)
tallies <- tallyVariants(bamfl, tally.Param, BPPARAM = bpp)
calling.filters <- VariantCallingFilters()
post.filters <- VariantPostFilters()
snp <- callVariants(tallies, calling.filters, post.filters)
sampleNames(snp) <- file_path_sans_ext(basename(bamfl))
mcols(snp) <- NULL
vcf <- asVCF(snp)
writeVcf(vcf, paste0(file_path_sans_ext(basename(bamfl)), ".vcf"), index = TRUE)
vcf <- readVcf(paste0(file_path_sans_ext(basename(bamfl)), ".vcf.bgz"))

##Due to constraints in memory, rsIDs are obtained on a chromosome by chromosome
##basis
##The resulting data frame is sorted as per the vcf and the rsIDs are injected
##into the vcf before writing it as a file again
all_snps <- SNPlocs.Hsapiens.dbSNP144.GRCh37
vcfID <- data.frame()
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
vcfID$chr <-  as.integer(gsub("chr", "", vcfID$chromosome)) 
dd <- vcfID[order(vcfID$chr, vcfID$Variant),]
dd$rsID[is.na(dd$rsID)] <- "."
names(vcf) <- dd$rsID

writeVcf(vcf, paste0(file_path_sans_ext(basename(bamfl)), ".vcf"), index = TRUE)