#' vizcnv
#'
#' @param cnvdata Path to the CNVs data file
#' @param species Species name from the following list: "Homo Sapiens", "Mus Musculus","Arabidopsis Thaliana" ,
#' "Drosophila Melanogaster", "Danio rerio", "Rattus norvegicus", or "Saccharomyces Cerevisiae" .
#'

vizcnv<- function( cnvdata, species) {
  # species
  if (species == "Homo Sapiens"){
    library(BSgenome.Hsapiens.UCSC.hg38)
    genome_species <- "hg38"
  } else if (species == "Mus Musculus"){
    library(BSgenome.Mmusculus.UCSC.mm10)
    genome_species <- "mm10"

  } else if (species == "Arabidopsis Thaliana"){
    library(BSgenome.Athaliana.TAIR.TAIR9)
    genome_species <- "TAIR9"

  } else if (species == "Caenorhabditis Elegans"){
    library(BSgenome.Celegans.UCSC.ce11)
    genome_species <- "ce2"

  } else if (species == "Drosophila Melanogaster"){
    library(BSgenome.Dmelanogaster.UCSC.dm6)
    genome_species <- "dm3"

  } else if (species == "Danio rerio"){
    library(BSgenome.Drerio.UCSC.danRer11)
    genome_species <- "danRer10"

  } else if (species == "Rattus_norvegicus"){
    library(BSgenome.Rnorvegicus.UCSC.rn5)
    genome_species <- "rn5"

  } else if (species == "Saccharomyces Cerevisiae"){
    library(BSgenome.Scerevisiae.UCSC.sacCer3)
    genome_species <- "sacCer3"

  } #else {
   # print("Non-valid species!")
  #}

  # load requirements
  library(shinyWidgets)
  library(shiny)
  library(readr)
  library(shinythemes)
  library(CNVRanger)
  library(Gviz)
  library(AnnotationHub)
  library(regioneR)
  library(shinycssloaders)
  library(dplyr)
  #read data
  calls <- read.csv(cnvdata, as.is=TRUE)
  #run app

  ui <- bootstrapPage(

    navbarPage(theme = shinytheme("flatly"), collapsible = TRUE,
               HTML('<a style="text-decoration:none;cursor:default;color:#FFFFFF;" class="active" href="#">CNVviz</a>'), id="nav",
               windowTitle = "CNVviz",
               tabPanel("Recurrent region",

                        sidebarLayout(
                          sidebarPanel(
                            pickerInput("chrI", span(shiny::tags$i(h6("Define the chromosome")), style="color:#045a8d"),
                                        if (species == "Homo Sapiens"){
                                          choices = c("chr1", "chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                                      "chr11", "chr12","chr13","chr14","chr15","chr16","chr17","chr18",
                                                      "chr19", "chr20", "chr21", "chr22", "chr23", "chrY", "chrX")
                                        } else if (species == "Mus Musculus"){
                                          choices = c("chr1", "chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                                      "chr11", "chr12","chr13","chr14","chr15","chr16","chr17","chr18",
                                                      "chr19", "chrY", "chrX")
                                        } else if (species == "Arabidopsis Thaliana"){
                                          choices = c("chr1", "chr2","chr3","chr4","chr5")
                                        }
                                        else if (species == "Drosophila Melanogaster"){
                                          choices = c("chr1", "chr2","chr3","chr4", "chrY", "chrX")

                                        } else if (species == "Danio rerio"){
                                          choices = c("chr1", "chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                                      "chr11", "chr12","chr13","chr14","chr15","chr16","chr17","chr18",
                                                      "chr19", "chr23", "chr24", "chr25")
                                        } else if (species == "Rattus_norvegicus"){
                                          choices = c("chr1", "chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                                      "chr11", "chr12","chr13","chr14","chr15","chr16","chr17","chr18",
                                                      "chr19", "chr20", "chrY", "chrX")
                                        } else if (species == "Saccharomyces Cerevisiae"){
                                          choices = c("chr1", "chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                                      "chr11", "chr12","chr13","chr14","chr15","chr16")
                                        } else {
                                          print("warning")
                                        },
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

                                        if (species == "Homo Sapiens"){
                                          choices = c("1", "2","3","4","5","6","7","8","9","10",
                                                      "11", "12","13","14","15","16","17","18",
                                                      "19", "20", "21", "22", "23", "Y", "X")
                                        } else if (species == "Mus Musculus"){
                                          choices = c("1", "2","3","4","5","6","7","8","9","10",
                                                      "11", "12","13","14","15","16","17","18",
                                                      "19", "Y", "X")
                                        } else if (species == "Arabidopsis Thaliana"){
                                          choices = c("1", "2","3","4","5")
                                        }
                                        else if (species == "Drosophila Melanogaster"){
                                          choices = c("1", "2","3","4", "Y", "X")
                                        } else if (species == "Danio rerio"){
                                          choices = c("1", "2","3","4","5","6","7","8","9","10",
                                                      "11", "12","13","14","15","16","17","18",
                                                      "19", "23", "24", "25")
                                        } else if (species == "Rattus_norvegicus"){
                                          choices = c("1", "2","3","4","5","6","7","8","9","10",
                                                      "11", "12","13","14","15","16","17","18",
                                                      "19", "20", "Y", "X")
                                        } else if (species == "Saccharomyces Cerevisiae"){
                                          choices = c("1", "2","3","4","5","6","7","8","9","10",
                                                      "11", "12","13","14","15","16")
                                        } else {
                                          print("warning")
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

                                        if (species == "Homo Sapiens"){
                                          choices = c("1", "2","3","4","5","6","7","8","9","10",
                                                      "11", "12","13","14","15","16","17","18",
                                                      "19", "20", "21", "22", "23", "Y", "X")
                                        } else if (species == "Mus Musculus"){
                                          choices = c("1", "2","3","4","5","6","7","8","9","10",
                                                      "11", "12","13","14","15","16","17","18",
                                                      "19", "Y", "X")
                                        } else if (species == "Arabidopsis Thaliana"){
                                          choices = c("1", "2","3","4","5")
                                        }
                                        else if (species == "Drosophila Melanogaster"){
                                          choices = c("1", "2","3","4", "Y", "X")
                                        } else if (species == "Danio rerio"){
                                          choices = c("1", "2","3","4","5","6","7","8","9","10",
                                                      "11", "12","13","14","15","16","17","18",
                                                      "19", "23", "24", "25")
                                        } else if (species == "Rattus_norvegicus"){
                                          choices = c("1", "2","3","4","5","6","7","8","9","10",
                                                      "11", "12","13","14","15","16","17","18",
                                                      "19", "20", "Y", "X")
                                        } else if (species == "Saccharomyces Cerevisiae"){
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
    ))

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
    ahDb <- AnnotationHub::query(ah, pattern = c(species, "EnsDb"))
    if (species == "Homo Sapiens"){
      ahEdb <- ahDb[["AH98047"]]
    } else if (species == "Mus Musculus"){
      ahEdb <- ahDb[[""]]
    } else if (species == "Arabidopsis Thaliana"){
      ahEdb <- ahDb[[""]]
    }
    else if (species == "Drosophila Melanogaster"){
      ahEdb <- ahDb[[""]]
    } else if (species == "Danio rerio"){
      ahEdb <- ahDb[[""]]
    } else if (species == "Rattus_norvegicus"){
      ahEdb <- ahDb[[""]]
    } else if (species == "Saccharomyces Cerevisiae"){
      ahEdb <- ahDb[[""]]
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
        plotRecurrentRegions(cnvrs, genome =genome_species, chr = input$chrI )   ###specie
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
