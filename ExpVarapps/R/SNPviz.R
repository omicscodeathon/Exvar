#' SNPviz
#'
#' @param datafile path to the SNP data vcf file
#'
#' @return
#' @export
#'
#' @examples
SNPviz <- function(datafile){

  library(shinyWidgets)
  library(shiny)
  library(readr)
  library(shinythemes)
  library(DT)
  library(echarts4r)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(GenomicFeatures)
  library(stringr)
  library(VariantAnnotation)
  library(data.table)
  library(VariantAnnotation)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(SNPlocs.Hsapiens.dbSNP144.GRCh37)

  snp_data <- fread(file=datafile, sep='\t', header=TRUE, skip='#CHROM')
  vcf <- readVcf(datafile, "hg19")
  vcf2 <- readVcf(datafile)
  ui <- bootstrapPage(

    navbarPage(theme = shinytheme("flatly"), collapsible = TRUE,
               HTML('<a style="text-decoration:none;cursor:default;color:#FFFFFF;" class="active" href="#">SNPviz</a>'), id="nav",
               windowTitle = "SNPviz",
               tabPanel("SNPs distribution",
                        box("SNPs types", downloadButton("downloadrplot", "Download the plot"),
                            echarts4rOutput("pie_snp_type", height = 500)),
                        box("TiTv distribution.", downloadButton("downloadrplot", "Download the plot"),
                            br(),
                            plotOutput("TiTv", height = 500))
               ),
               tabPanel("Trinucleotide Pattern",
                        box("Trinucleotide pattern plot. ", downloadButton("downloadrplot", "Download the plot"),
                            plotOutput("Tri_pattern", height = 500))
               ),
               tabPanel("Amino Acid Changes",
                        box("Trinucleotide pattern plot.", downloadButton("downloadrplot", "Download the plot"),
                            echarts4rOutput("Variantions", height = 500))

               )
    ))


  # Define server logic required to draw a histogram
  server <- function(input, output) {

    # detect snps
    snp_data$snp <- "Normal"
    snp_data$snp[snp_data$REF == "A" & snp_data$ALT == "T"] <- "A>T"
    AtoT <-nrow(snp_data[snp_data$snp == "A>T",])
    snp_data$snp[snp_data$REF == "A" & snp_data$ALT == "C"] <- "A>C"
    AtoC <-nrow(snp_data[snp_data$snp == "A>C",])
    snp_data$snp[snp_data$REF == "A" & snp_data$ALT == "G"] <- "A>G"
    AtoG <-nrow(snp_data[snp_data$snp == "A>G",])
    snp_data$snp[snp_data$REF == "T" & snp_data$ALT == "A"] <- "T>A"
    TtoA <-nrow(snp_data[snp_data$snp == "T>A",])
    snp_data$snp[snp_data$REF == "T" & snp_data$ALT == "C"] <- "T>C"
    TtoC <-nrow(snp_data[snp_data$snp == "T>C",])
    snp_data$snp[snp_data$REF == "T" & snp_data$ALT == "G"] <- "T>G"
    TtoG <-nrow(snp_data[snp_data$snp == "T>G",])
    snp_data$snp[snp_data$REF == "C" & snp_data$ALT == "T"] <- "C>T"
    CtoT <-nrow(snp_data[snp_data$snp == "C>T",])
    snp_data$snp[snp_data$REF == "C" & snp_data$ALT == "A"] <- "C>A"
    CtoA <-nrow(snp_data[snp_data$snp == "C>A",])
    snp_data$snp[snp_data$REF == "C" & snp_data$ALT == "G"] <- "C>G"
    CtoG <-nrow(snp_data[snp_data$snp == "C>G",])
    snp_data$snp[snp_data$REF == "G" & snp_data$ALT == "T"] <- "G>T"
    GtoT <-nrow(snp_data[snp_data$snp == "G>T",])
    snp_data$snp[snp_data$REF == "G" & snp_data$ALT == "C"] <- "G>C"
    GtoC<-nrow(snp_data[snp_data$snp == "G>C",])
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

    #################
    rd <- rowRanges(vcf)
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
    ## plot
    output$Tri_pattern <- renderPlot({
      ggplot(tbl_dat,aes(x=Var1,y=Freq,fill=Var1))+
        geom_bar(stat='identity')+
        labs(x="",y="Variants",fill="")+
        theme(legend.position = "none")

    })

    #######
    # filter a chromosome  ---------------------------
    chr = "chr1"
    vcf_chr1 <- vcf[grepl(names(vcf),pattern = chr)]
    rd_chr1 <- rowRanges(vcf_chr1)
    # Retrive dbSNP data
    all_snps <- SNPlocs.Hsapiens.dbSNP144.GRCh37
    all_snps
    tar_chr <- as.vector(seqnames(rd_chr1)@values)
    tar_chr <- gsub("chr", "", tar_chr)
    tar_chr[grepl(tar_chr, pattern = "M")] <- "MT"
    my_snps <- snpsBySeqname(all_snps, c(tar_chr))
    my_snps[1:2]
    #Convert seqInfo to UCSC style
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
    #####
    chr <- 1:22
    chrxym <- c("X", "Y", "M")
    chr <- append(chr, chrxym)
    mainChromosomes <- paste0("chr", chr)
    seqlevels(vcf2) <- mainChromosomes
    vcf2[grepl("MT", vcf)] <- "M"
    genome(vcf2) <- "hg19"
    seqlengths(vcf2) <- head(seqlengths(TxDb.Hsapiens.UCSC.hg19.knownGene), 25L)
    vcf_coding <- predictCoding(vcf2, TxDb.Hsapiens.UCSC.hg19.knownGene,
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

  }

  # Run the application
  shinyApp(ui = ui, server = server)

}
