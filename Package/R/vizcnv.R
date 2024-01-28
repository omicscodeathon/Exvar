#' vizcnv
#'
#' @param cnvdata Path to the CNVs data file
#' @export
#'

vizcnv<- function( cnvdata) {
# specie
cat(paste0("These are the species currently supported by ExpVar: \n",
           "[1] Homo sapiens \n",
           "[2] Mus musculus \n",
           "[3] Drosophila melanogaster \n",
           "[4] Danio rerio \n",
           "[5] Rattus norvegicus \n",
           "[6] Saccharomyces cerevisiae \n"))
species <- readline("Type the number of the species that you would like to use as a reference: ")
#  load requirements
library(shiny)
library(shinyWidgets)
library(shinythemes)
library(shinycssloaders)
library(readr)
library(CNVRanger)
library(Gviz)
library(AnnotationHub)
library(regioneR)
library(dplyr)
# species specific requirements
switch(species,
       "1"={
         library(BSgenome.Hsapiens.UCSC.hg38)
         genome_species <- "hg38"
         sp <- "Homo sapiens"
       },
       "2"={
         library(BSgenome.Mmusculus.UCSC.mm10)
         genome_species <- "mm10"
         sp <- "Mus musculus"
       },
       "3"={
         library(BSgenome.Dmelanogaster.UCSC.dm6)
         genome_species <- "dm3"
         sp <-"Drosophila melanogaster"
       },
       "4"={
         library(BSgenome.Drerio.UCSC.danRer11)
         genome_species <- "danRer10"
         sp <- "Danio rerio"
       },
       "5"={ 
         library(BSgenome.Rnorvegicus.UCSC.rn5)
         genome_species <- "rn5"
         sp <- "Rattus norvegicus"
       },
       "6"={
         library(BSgenome.Scerevisiae.UCSC.sacCer3)
         genome_species <- "sacCer3"
         sp <-"Saccharomyces cerevisiae"
       }
)

  #read data
  calls <- read.csv(cnvdata, as.is=TRUE)
  #run app
  
  ui <- navbarPage(
    theme = shinytheme("flatly"),
    HTML('<a style="text-decoration:none;
               cursor:default;
                    color:white;
                    " class="active" href="#">Visualize CNV data</a>'),
    windowTitle ="vizcnv",
               tabPanel("Recurrent region",
                        
                        sidebarLayout(
                          sidebarPanel(
                            pickerInput("chrI", span(shiny::tags$i(h6("Define the chromosome")), style="color:#045a8d"),
                                        if (species == "1"){
                                          choices = c("chr1", "chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                                      "chr11", "chr12","chr13","chr14","chr15","chr16","chr17","chr18",
                                                      "chr19", "chr20", "chr21", "chr22", "chr23", "chrY", "chrX")
                                        } else if (species == "2"){
                                          choices = c("chr1", "chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                                      "chr11", "chr12","chr13","chr14","chr15","chr16","chr17","chr18",
                                                      "chr19", "chrY", "chrX")
                                        }
                                        else if (species == "3"){
                                          choices = c("chr1", "chr2","chr3","chr4", "chrY", "chrX")
                                          
                                        } else if (species == "4"){
                                          choices = c("chr1", "chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                                      "chr11", "chr12","chr13","chr14","chr15","chr16","chr17","chr18",
                                                      "chr19", "chr23", "chr24", "chr25")
                                        } else if (species == "5"){
                                          choices = c("chr1", "chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                                      "chr11", "chr12","chr13","chr14","chr15","chr16","chr17","chr18",
                                                      "chr19", "chr20", "chrY", "chrX")
                                        } else if (species == "6"){
                                          choices = c("chr1", "chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                                      "chr11", "chr12","chr13","chr14","chr15","chr16")
                                        } ,
                                        selected = c("chr1"),
                                        multiple = FALSE),
                            br(),
                            
                            sliderInput("pvalue",
                                        span(shiny::tags$i(h6("Define the significance threshold to filter recurrent CNVs.")), style="color:#045a8d"),
                                        min = 0,
                                        max = 1,
                                        value=0.05 ),
                            
                            br(),
                            downloadButton("downloadrplot", "Download the plot")
                          ),
                          mainPanel( shinycssloaders::withSpinner(plotOutput("recurrent_plot"))
                          )
                        )),
               
               tabPanel("Overlap analysis",
                        sidebarLayout(
                          sidebarPanel(
                            pickerInput("chrII", span(shiny::tags$i(h6("Define the chromosome")), style="color:#045a8d"),
                                        
                                        if (species == "1"){
                                          choices = c("1", "2","3","4","5","6","7","8","9","10",
                                                      "11", "12","13","14","15","16","17","18",
                                                      "19", "20", "21", "22", "23", "Y", "X")
                                        } else if (species == "2"){
                                          choices = c("1", "2","3","4","5","6","7","8","9","10",
                                                      "11", "12","13","14","15","16","17","18",
                                                      "19", "Y", "X")
                                        } else if (species == "3"){
                                          choices = c("1", "2","3","4", "Y", "X")
                                        } else if (species == "4"){
                                          choices = c("1", "2","3","4","5","6","7","8","9","10",
                                                      "11", "12","13","14","15","16","17","18",
                                                      "19", "23", "24", "25")
                                        } else if (species == "5"){
                                          choices = c("1", "2","3","4","5","6","7","8","9","10",
                                                      "11", "12","13","14","15","16","17","18",
                                                      "19", "20", "Y", "X")
                                        } else if (species == "6"){
                                          choices = c("1", "2","3","4","5","6","7","8","9","10",
                                                      "11", "12","13","14","15","16")
                                        },
                                        selected = c("1"),
                                        multiple = FALSE),
                            br(),
                            pickerInput("biotype",  span(shiny::tags$i(h6("Define the target gene biotype.")), style="color:#045a8d"),
                                        choices = c("protein_coding", "pseudogene", "rRNA", "miRNA","misc_RNA","snRNA", "snoRNA"),
                                        selected = c("protein_coding"),
                                        multiple = FALSE),
                            br(),
                            downloadButton("downloadoplot", "Download  the plot")
                            
                            
                          ),
                          mainPanel( shinycssloaders::withSpinner( plotOutput("overlap_plot")))
                          
                        )),
               
               tabPanel("Overlap permutation test",
                        sidebarLayout(
                          sidebarPanel(
                            pickerInput("chrIII", span(shiny::tags$i(h6("Define the chromosome")), style="color:#045a8d"),
                                        
                                        if (species == "1"){
                                          choices = c("1", "2","3","4","5","6","7","8","9","10",
                                                      "11", "12","13","14","15","16","17","18",
                                                      "19", "20", "21", "22", "23", "Y", "X")
                                        } else if (species == "2"){
                                          choices = c("1", "2","3","4","5","6","7","8","9","10",
                                                      "11", "12","13","14","15","16","17","18",
                                                      "19", "Y", "X")
                                        } else if (species == "3"){
                                          choices = c("1", "2","3","4", "Y", "X")
                                        } else if (species == "4"){
                                          choices = c("1", "2","3","4","5","6","7","8","9","10",
                                                      "11", "12","13","14","15","16","17","18",
                                                      "19", "23", "24", "25")
                                        } else if (species == "5"){
                                          choices = c("1", "2","3","4","5","6","7","8","9","10",
                                                      "11", "12","13","14","15","16","17","18",
                                                      "19", "20", "Y", "X")
                                        } else if (species == "6"){
                                          choices = c("1", "2","3","4","5","6","7","8","9","10",
                                                      "11", "12","13","14","15","16")
                                        } else {
                                          print("warning")
                                        },
                                        selected = c("1"),
                                        multiple = FALSE),
                            br(),
                            pickerInput("biotypeI",  span(shiny::tags$i(h6("Define the target gene biotype.")), style="color:#045a8d"),
                                        choices = c("protein_coding", "pseudogene", "rRNA", "miRNA","misc_RNA","snRNA", "snoRNA"),
                                        selected = c("protein_coding"),
                                        multiple = FALSE),
                            br(),
                            downloadButton("downloadpplot", "Download the plot")
                          ),
                          mainPanel( shinycssloaders::withSpinner( plotOutput("permutation_plot")))
                        ))
    )
  
  server <- function(input, output) {
    
    length(unique(calls[,"sample_ID"]))
    #Representation as a GRangesList
    grl <- makeGRangesListFromDataFrame(calls,
                                        split.field = "sample_ID", keep.extra.columns = TRUE)
    grl <- sort(grl)
    #Representation as a RaggedExperiment
    ra <- RaggedExperiment(grl)
    weight <- rnorm(ncol(ra), mean=1100, sd=100)
    fcr <- rnorm(ncol(ra), mean=6, sd=1.5)
    colData(ra)$weight <- round(weight, digits=2)
    colData(ra)$fcr <- round(fcr, digits=2)
    #Identifying recurrent regions
    cnvrs <- populationRanges(grl, density = 0.1, est.recur = TRUE)
    
    # plot1
    output$recurrent_plot <- renderPlot({
      subset(cnvrs, pvalue < input$pvalue)
      plotRecurrentRegions(cnvrs, genome = genome_species, chr = input$chrI )
    })
    
    #Overlap analysis of CNVs with functional genomic regions
    ah <- AnnotationHub::AnnotationHub()
    ahDb <- AnnotationHub::query(ah, pattern = c(sp, "EnsDb"))
    ahDb <- AnnotationHub::query(ah, pattern = c(sp))
    if (species == "1"){
      ahEdb <- ahDb[["AH104864"]]
    } else if (species == "2"){
      ahEdb <- ahDb[["AH104895"]]
    } else if (species == "3"){
      ahEdb <- ahDb[["AH104833"]]
    } else if (species == "4"){
      ahEdb <- ahDb[["AH104837"]]
    } else if (species == "5"){
      ahEdb <- ahDb[["AH104966"]]
    } else if (species == "6"){
      ahEdb <- ahDb[["AH104972"]]
    }
    bt.genes <- ensembldb::genes(ahEdb)
    GenomeInfoDb::seqlevelsStyle(bt.genes) <- "UCSC"
    
    # plot 2
    output$overlap_plot <- renderPlot({
      sel.genes <- subset(bt.genes, seqnames %in% paste0("chr", input$chrII : input$chrII))
      sel.genes <- subset(sel.genes, gene_biotype == input$biotype)
      sel.cnvrs <- subset(cnvrs, seqnames %in% paste0("chr", input$chrII : input$chrII))
      #Finding and illustrating overlaps
      olaps <- GenomicRanges::findOverlaps(sel.genes, sel.cnvrs, ignore.strand=TRUE)
      qh <- S4Vectors::queryHits(olaps)
      sh <- S4Vectors::subjectHits(olaps)
      cgenes <- sel.genes[qh]
      cgenes$type <- sel.cnvrs$type[sh]
      subset(cgenes, select = "type")
      #plot
      cnvOncoPrint(grl, cgenes)
    })
    
    #Overlap permutation test
    output$permutation_plot <- renderPlot({
      sel.genes <- subset(bt.genes, seqnames %in% paste0("chr", input$chrIII : input$chrIII))
      sel.genes <- subset(sel.genes, gene_biotype == input$biotypeI)
      sel.cnvrs <- subset(cnvrs, seqnames %in% paste0("chr", input$chrIII : input$chrIII))
      
      res <- suppressWarnings(overlapPermTest(A=sel.cnvrs, B=sel.genes, ntimes=100,
                                              genome=genome_species, mask=NA, per.chromosome=TRUE, count.once=TRUE))
      #plot
      plot(res)
    })
    ####download plots
    #recurrent pot
    output$downloadrplot <- downloadHandler(
      filename = function(){
        paste("recurrent_plot", "png", sep=".")
      },
      content = function(file){
        png(file)
        plotRecurrentRegions(cnvrs, genome =genome_species, chr = input$chrI )   
        dev.off()
      })
    #
    output$downloadoplot <- downloadHandler(
      filename = function(){
        paste("Overlap_analysis_plot", "png", sep=".")
      },
      content = function(file){
        png(file)
        grl <- makeGRangesListFromDataFrame(calls,
                                            split.field = "sample_ID", keep.extra.columns = TRUE)
        grl <- sort(grl)
        sel.genes <- subset(bt.genes, seqnames %in% paste0("chr", input$chrII : input$chrII))
        sel.genes <- subset(sel.genes, gene_biotype == input$biotype)
        sel.cnvrs <- subset(cnvrs, seqnames %in% paste0("chr", input$chrII : input$chrII))
        olaps <- GenomicRanges::findOverlaps(sel.genes, sel.cnvrs, ignore.strand=TRUE)
        qh <- S4Vectors::queryHits(olaps)
        sh <- S4Vectors::subjectHits(olaps)
        cgenes <- sel.genes[qh]
        cgenes$type <- sel.cnvrs$type[sh]
        subset(cgenes, select = "type")
        cnvOncoPrint(grl, cgenes)
        dev.off()
      })
    #
    output$downloadpplot <- downloadHandler(
      filename = function(){
        paste("Overlap_permutation_test_plot", "png", sep=".")
      },
      content = function(file){
        png(file)
        sel.genes <- subset(bt.genes, seqnames %in% paste0("chr", input$chrIII : input$chrIII))
        sel.genes <- subset(sel.genes, gene_biotype == input$biotypeI)
        sel.cnvrs <- subset(cnvrs, seqnames %in% paste0("chr", input$chrIII : input$chrIII))
        
        res <- suppressWarnings(overlapPermTest(A=sel.cnvrs, B=sel.genes, ntimes=100,
                                                genome=genome_species, mask=NA, per.chromosome=TRUE, count.once=TRUE))
        
        plot(res)
        dev.off()
      })
  }
  
  # Run the application
  shinyApp(ui = ui, server = server)

}
