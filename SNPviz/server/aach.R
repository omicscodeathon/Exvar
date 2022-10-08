vcf <- readVcf("demo/BT_snp.vcf")

#
chr <- 1:22
chrxym <- c("X", "Y", "M")
chr <- append(chr, chrxym)
mainChromosomes <- paste0("chr", chr)
seqlevels(vcf) <- mainChromosomes
vcf[grepl("MT", vcf)] <- "M"
genome(vcf) <- "hg19"
seqlengths(vcf) <- head(seqlengths(TxDb.Hsapiens.UCSC.hg19.knownGene), 25L)
vcf_coding <- predictCoding(vcf, TxDb.Hsapiens.UCSC.hg19.knownGene, 
                            BSgenome.Hsapiens.UCSC.hg19)
#####
coding <- vcf_coding

coding[1]
matA <- data.frame(Variant=names(coding),
                   chromosome=seqnames(coding),
                   start=start(coding),end=end(coding),
                   ref_allele=as.character(coding$REF),
                   alt_allele=unlist(lapply(lapply(
                     coding$ALT,`[[`,1),as.character)),
                   GeneID=coding$GENEID,
                   TxID=coding$TXID,
                   Protein_posi=unlist(lapply(lapply(
                     coding$PROTEINLOC,`[[`,1),as.integer)),
                   ref_AA=as.character(coding$REFAA),
                   alt_AA=as.character(coding$VARAA),
                   Type=coding$CONSEQUENCE)


matA$aaChange <- paste0("p.",matA$ref_AA,matA$Protein_posi,matA$alt_AA)
matA <- dplyr::select(matA,-Protein_posi,-ref_AA,-alt_AA)

#Annotation table ~ Amino Acid Changes
matA[1:2, ]

#How many variations in coding region
var_in_coding <- data.frame(varName = names(vcf_chr1), in_coding = names(vcf_chr1) %in% 
                              matA$Variant, stringsAsFactors = FALSE)
table(var_in_coding$in_coding)
#How many types of mutations in coding region
taC <- table(matA$Type)
taC_dat <- as.data.frame(taC)
taC

##
output$aach<- renderPlot({
  ggplot(taC_dat,aes(x=Var1,y=Freq,fill=Var1))+
    geom_bar(stat='Identity')+
    labs(x="",y="Counts",fill="")+
    theme(legend.position = "none")
})


#Integrate SNP and amino acid change into single table
matS$GeneID <- matA$GeneID[match(matS$Variant, matA$Variant)]
matS$AAChange <- matA$GeneID[match(matS$Variant, matA$Variant)]
matS[1:2, ]
