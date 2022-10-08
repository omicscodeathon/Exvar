#visualize the vcf file
#library(VariantAnnotation)
#library(BSgenome.Hsapiens.UCSC.hg19)
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#library(SNPlocs.Hsapiens.dbSNP144.GRCh37)

vcf <- readVcf("demo/all_variants.vcf", "hg19")
vcf  # output

############## ---------------------------      genotype types         ---------------------------#############

# genotype types
tbl <- table(geno(vcf)$GT)
tbl_dat <- as.data.frame(tbl)
tbl

# gt types plot 
ggplot(tbl_dat,aes(x=Var1,y=Freq,fill=Var1))+
  geom_bar(stat='Identity')+
  labs(x="",y="Counts",fill="")+
  theme_classic()


tmp_vcf_data = df
#vcf = collaps

## Expand the CollapsedVCF to ExpandedVCF
evcf <- expand(vcf)
transition <- isTransition(evcf)
transition
DataFrame(ref(evcf)[transition], alt(evcf)[transition])

## S4 method for signature 'ExpandedVCF'============== filter SNPs
## Expand the CollapsedVCF to ExpandedVCF
evcf <- expand(vcf)
SNV = isSNV(evcf)
DataFrame(ref(evcf)[SNV], alt(evcf)[SNV])

## S4 method for signature 'ExpandedVCF'
isInsertion(x, ...)
## S4 method for signature 'ExpandedVCF'
isDeletion(x, ...)
## S4 method for signature 'ExpandedVCF' ============== filter indels
indel = isIndel(evcf)
DataFrame(ref(evcf)[indel], alt(evcf)[indel])

## S4 method for signature 'ExpandedVCF'
isDelins(x, ...)
## S4 method for signature 'ExpandedVCF'
isTransition(x, ...)
## S4 method for signature 'ExpandedVCF'
isSubstitution(x, ...)

##
isInsertion(x, ..., singleAltOnly = TRUE)
isSNV(df, singleAltOnly=FALSE)


# filte cnv ????????


############## ---------------------------       detect variant type       ---------------------------#############
# differentiate SNP/INS/DEL/Others
for (k in 1:length(tmp_vcf_data$variant)) {
  if (width(tmp_vcf_data$REF[k]) < width(tmp_vcf_data$ALT[k])) {
    tmp_vcf_data$mutType[k] <- "INS"
  } else if (width(tmp_vcf_data$REF[k]) > width(tmp_vcf_data$ALT[k])) {
    tmp_vcf_data$mutType[k] <- "DEL"
  } else if (width(tmp_vcf_data$REF[k]) == 1 & width(tmp_vcf_data$ALT[k]) == 1) {
    tmp_vcf_data$mutType[k] <- "SNP"
  } else {
    tmp_vcf_data$mutType[k] <- "Others"
  }
}
# 
tbl <- table(tmp_vcf_data$mutType)
tbl_dat <- as.data.frame(tbl)
tbl

# plot variants types
library(ggplot2)
ggplot(tbl_dat,aes(x=Var1,y=Freq,fill=Var1))+
  geom_bar(stat = 'identity')+
  labs(x="",y="Mutations",fill="")+
  theme_classic()
############## ---------------------------      + keep only snps        ---------------------------#############











