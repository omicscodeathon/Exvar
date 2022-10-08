
library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37) 

vcf <- readVcf("demo/BT_snp.vcf", "hg19")

# filter a chromosome  ---------------------------
chr = "chr1"   # input

# create the loop for all chromosoms
#if (chr = "chr1" ){
#  chr2 <- "chr1:"
#} else if (chr="chr2"){
#  chr2 <- "chr2:"
#}
#
vcf_chr1 <- vcf[grepl(names(vcf),pattern = chr)]
rd_chr1 <- rowRanges(vcf_chr1)

# Retrive dbSNP data
all_snps <- SNPlocs.Hsapiens.dbSNP144.GRCh37
all_snps

# ------------
#Retrieve SNPs by chromosome

#tar_chr <- as.vector(seqnames(rd_chr1)@values)
#my_snps <- snpsBySeqname(all_snps, c(tar_chr))  ## eror

#genome in rd_chr1 is hg19, but in all_snps is GRCh37.p13
#Unify seqnames and Process
tar_chr <- as.vector(seqnames(rd_chr1)@values)
tar_chr <- gsub("chr", "", tar_chr)
tar_chr[grepl(tar_chr, pattern = "M")] <- "MT"
my_snps <- snpsBySeqname(all_snps, c(tar_chr))
my_snps[1:2]

#Convert seqInfo to UCSC style   # performed only if needed
## change seqlevelsStyle
seqlevelsStyle(my_snps) <- "UCSC"
# change genome
genome(my_snps) <- "hg19"

# make rsID table
snp_ID <- data.frame(posIDX = paste0(seqnames(my_snps), ":", pos(my_snps)), rsID = my_snps$RefSNP_id, stringsAsFactors = FALSE)
head(snp_ID)

#Generate Variant table
matV1 <- data.frame(Variant = names(rd_chr1), stringsAsFactors = FALSE)
matV1[1:2, ]
#Extract information of variants
matV1$chromosome <- gsub("(.*):(.*)_(.*)/(.*)", "\\1", matV1$Variant)
matV1$start <- gsub("(.*):(.*)_(.*)/(.*)", "\\2", matV1$Variant)
matV1$end <- gsub("(.*):(.*)_(.*)/(.*)", "\\2", matV1$Variant)
matV1$ref_allele <- gsub("(.*):(.*)_(.*)/(.*)", "\\3", matV1$Variant)
matV1$alt_allele <- gsub("(.*):(.*)_(.*)/(.*)", "\\4", matV1$Variant)
matV1$posIDX <- gsub("(.*)_(.*)", "\\1", matV1$Variant)
matV1[1:2, ]

#Annotation table ~ SNP_ID
matS <- merge(matV1,snp_ID,all.x=TRUE,by="posIDX")
matS <- dplyr::select(matS,-posIDX)
matS[1:2,]
# How many variations in dbSNP
taC2 <- table(!is.na(matS$rsID))
taC2_dat <- as.data.frame(taC2)
taC2


#Variations in dbSNP ~ Plotting

output$Variantions <- renderPlot({
  ggplot(taC2_dat,aes(x=Var1,y=Freq,fill=Var1))+
    geom_bar(stat='Identity')+
    labs(x="",y="Counts",fill="in_dbSNP")+
    theme(legend.position = "none")
  
  
})





















