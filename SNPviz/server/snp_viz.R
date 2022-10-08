
library(echarts4r)
library(readr)
library(dplyr)
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicFeatures)
library(stringr)
library(VariantAnnotation)
library(data.table)
df <- fread(file='demo/BT_snp.vcf', sep='\t', header=TRUE, skip='#CHROM')
snp_data= df
#---------------------------------------------------------------------------------------------------

# filtering
## filter a chromosome
#chr = "chr2"
#data<- subset(snp_data, CHROM == chr) 


############### ------------------------  PART1 ------------------------  ####################
# detect snps
snp_data$snp <- "Normal"
# A>T
snp_data$snp[snp_data$REF == "A" & snp_data$ALT == "T"] <- "A>T"
AtoT <-nrow(snp_data[snp_data$snp == "A>T",])
# A>C
snp_data$snp[snp_data$REF == "A" & snp_data$ALT == "C"] <- "A>C"
AtoC <-nrow(snp_data[snp_data$snp == "A>C",])
# A>G
snp_data$snp[snp_data$REF == "A" & snp_data$ALT == "G"] <- "A>G"
AtoG <-nrow(snp_data[snp_data$snp == "A>G",])

# T<A
snp_data$snp[snp_data$REF == "T" & snp_data$ALT == "A"] <- "T>A"
TtoA <-nrow(snp_data[snp_data$snp == "T>A",])
# T>C
snp_data$snp[snp_data$REF == "T" & snp_data$ALT == "C"] <- "T>C"
TtoC <-nrow(snp_data[snp_data$snp == "T>C",])
# T>G
snp_data$snp[snp_data$REF == "T" & snp_data$ALT == "G"] <- "T>G"
TtoG <-nrow(snp_data[snp_data$snp == "T>G",])

# c>T
snp_data$snp[snp_data$REF == "C" & snp_data$ALT == "T"] <- "C>T"
CtoT <-nrow(snp_data[snp_data$snp == "C>T",])
# c>A
snp_data$snp[snp_data$REF == "C" & snp_data$ALT == "A"] <- "C>A"
CtoA <-nrow(snp_data[snp_data$snp == "C>A",])
# C>G
snp_data$snp[snp_data$REF == "C" & snp_data$ALT == "G"] <- "C>G"
CtoG <-nrow(snp_data[snp_data$snp == "C>G",])

# G>T
snp_data$snp[snp_data$REF == "G" & snp_data$ALT == "T"] <- "G>T"
GtoT <-nrow(snp_data[snp_data$snp == "G>T",])
# G>C
snp_data$snp[snp_data$REF == "G" & snp_data$ALT == "C"] <- "G>C"
GtoC<-nrow(snp_data[snp_data$snp == "G>C",])
# G>A
snp_data$snp[snp_data$REF == "G" & snp_data$ALT == "A"] <- "G>A"
GtoA <-nrow(snp_data[snp_data$snp == "G>A",])

## creat dataframe
SNP_type <- c("AtoC", "AtoG", "AtoT", "CtoA", "CtoG", "CtoT","TtoC", "TtoG", "TtoA", "GtoA", "GtoC", "GtoT")
Value <- c(AtoC, AtoG, AtoT, CtoA, CtoG, CtoT, TtoC, TtoG, TtoA, GtoA, GtoC, GtoT)

SNPs_types_count <- data.frame(SNP_type, Value)

## Pie graph snp type ------------------------ 

output$pie_snp_type<-renderEcharts4r({
  SNPs_types_count |>   # change output command
    head() |> 
    e_charts(SNP_type) |> 
    e_pie(Value) |> 
    e_title("Pie chart") 
})

  ############### ------------------------  PART2 ------------------------  ####################
# identify Ti and Tv
snp_data2 = snp_data
# Transition (Ti) : purine-to-purine, pyrimidine-to-pyrimidine
ti <- c("A>G","G>A","C>T","T>C")
# Transveersion (Tv) : purine-to-pyrimidine, pyrimidine-to-purine
tv <- c("A>T","A>C","G>T","G>C","C>A","C>G","T>A","T>G")

snp_data2$TiTv[snp_data2$snp %in% ti] <- "Ti"
snp_data2$TiTv[snp_data2$snp %in% tv] <- "Tv"
snp_data2[1:2,]


# plot 
output$TiTv <- renderPlot({
  ggplot(as.data.frame(table(snp_data2$TiTv)),aes(x=Var1,y=Freq,fill=Var1))+
    geom_bar(stat = 'identity')+labs(x="",y="Mutations",fill="")+
    theme(legend.position = "none")
})

####
############### ------------------------  PART3 ------------------------  ####################
#### Trinucleotide motif analysis
vcf <- readVcf("demo/BT_snp.vcf", "hg19")
rd <- rowRanges(vcf)
#rd[1:2]
# 
rd_idx <- str_split(names(rd), "_", simplify = T)
rd_idx[1:2, ]

rd_sub <- rd[rd_idx[, 2] == "C/T"]
rd_sub[1:2, ]

# Extract sequences beneath the mutation from -1 to +1
rd_sub$triNu <- getSeq(Hsapiens,
                       seqnames(rd_sub),
                       start=start(rd_sub)-1,
                       end=end(rd_sub)+1)
rd_sub[1:2]

## Trinucleotide pattern

tbl <- table(rd_sub$triNu)
tbl_dat <- as.data.frame(tbl)
tbl

## plot
output$Tri_pattern <- renderPlot({
  ggplot(tbl_dat,aes(x=Var1,y=Freq,fill=Var1))+
    geom_bar(stat='identity')+
    labs(x="",y="Variants",fill="")+
    theme(legend.position = "none")
  
})








