#' Title
#'
#' @param datafile path to the CNVs csv file
#'
#' @return
#' @export
#'
#' @examples
CNVviz <- function(datafile) {

  library(shinyWidgets)
  library(shiny)
  library(readr)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(shinythemes)
  library(CNVRanger)
  library(Gviz)
  library(AnnotationHub)
  library(regioneR)

  inputdata<- read.csv(datafile, as.is=TRUE) = input
  ui <- bootstrapPage(

    navbarPage(theme = shinytheme("flatly"), collapsible = TRUE,
               HTML('<a style="text-decoration:none;cursor:default;color:#FFFFFF;" class="active" href="#">CNVviz</a>'), id="nav",
               windowTitle = "CNVviz",
               tabPanel("Recurrent region",

                        sidebarLayout(
                          sidebarPanel(

                            span(shiny::tags$i(h6("Please define the chromosome")), style="color:#045a8d"),

                            pickerInput("chrI", "Chromosome:",
                                        choices = c("chr1", "chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                                    "chr11", "chr12","chr13","chr14","chr15","chr16","chr17","chr18",
                                                    "chr19", "chr20", "chr21", "chr22", "chr23", "chrY", "chrX"),
                                        selected = c("chr1"),
                                        multiple = FALSE),
                            br(),
                            span(shiny::tags$i(h6("Please define the significance threshold to filter recurrent CNVs.")), style="color:#045a8d"),
                            br(),
                            sliderInput("pvalue",
                                        "P_value:",
                                        min = 0,
                                        max = 1,
                                        value=0.05,
                                        #timeFormat="%d %b"
                            ),

                            br(),
                            downloadButton("downloadrplot", "Download the plot")
                          ),
                          mainPanel(
                            plotOutput("recurrent_plot") )
                        )),

               tabPanel("Overlap analysis",
                        sidebarLayout(
                          sidebarPanel(

                            span(shiny::tags$i(h6("Define the chromosome")), style="color:#045a8d"),

                            pickerInput("chrII", "Chromosome:",
                                        choices = c("1", "2","3","4","5","6","7","8","9","10",
                                                    "11", "12","13","14","15","16","17","18",
                                                    "19", "20", "21", "22", "23", "Y", "X"),
                                        selected = c("1"),
                                        multiple = FALSE),
                            br(),
                            span(shiny::tags$i(h6("Define the target gene biotype.")), style="color:#045a8d"),
                            br(),
                            pickerInput("biotype", "Gene Biotype:",
                                        choices = c("protein_coding", "pseudogene", "rRNA", "miRNA","misc_RNA","snRNA", "snoRNA"),
                                        selected = c("protein_coding"),
                                        multiple = FALSE),
                            br(),
                            downloadButton("downloadoplot", "Download  the plot")


                          ),
                          mainPanel(plotOutput("overlap_plot") )
                        )),

               #
               tabPanel("Overlap permutation test",
                        sidebarLayout(
                          sidebarPanel(

                            span(shiny::tags$i(h6("Define the chromosome")), style="color:#045a8d"),

                            pickerInput("chrIII", "Chromosome:",
                                        choices = c("1", "2","3","4","5","6","7","8","9","10",
                                                    "11", "12","13","14","15","16","17","18",
                                                    "19", "20", "21", "22", "23", "Y", "X"),
                                        selected = c("1"),
                                        multiple = FALSE),
                            br(),
                            span(shiny::tags$i(h6("Define the target gene biotype.")), style="color:#045a8d"),
                            br(),
                            pickerInput("biotypeI", "Gene Biotype:",
                                        choices = c("protein_coding", "b2"),
                                        selected = c("protein_coding"),
                                        multiple = FALSE),
                            br(),
                            downloadButton("downloadpplot", "Download the plot")


                          ),
                          mainPanel(plotOutput("permutation_plot") )

                        ))


    ))


  # Define server logic required to draw a histogram
  server <- function(input, output) {
    # data

    length(unique(inputdata[,"NE_id"]))
    #Representation as a GRangesList
    grl <- makeGRangesListFromDataFrame(calls,
                                        split.field = "NE_id", keep.extra.columns = TRUE)
    grl <- sort(grl)

    #Representation as a RaggedExperiment
    ra <- RaggedExperiment(grl)
    assay(ra[1:5,1:5])

    weight <- rnorm(ncol(ra), mean=1100, sd=100)
    fcr <- rnorm(ncol(ra), mean=6, sd=1.5)
    colData(ra)$weight <- round(weight, digits=2)
    colData(ra)$fcr <- round(fcr, digits=2)
    #colData(ra)

    #Trimming low-density areas
    cnvrs <- populationRanges(grl, density = 0.1)
    #Identifying recurrent regions
    cnvrs <- populationRanges(grl, density = 0.1, est.recur = TRUE)

    # plot1
    output$recurrent_plot <- renderPlot({
      subset(cnvrs, pvalue < input$pvalue)
      plotRecurrentRegions(cnvrs, genome = "hg38", chr = input$chrI )
    })

    #Overlap analysis of CNVs with functional genomic regions

    ah <- AnnotationHub::AnnotationHub()
    ahDb <- AnnotationHub::query(ah, pattern = c("Homo sapiens", "EnsDb"))
    ahEdb <- ahDb[["AH53211"]]
    bt.genes <- ensembldb::genes(ahEdb)
    GenomeInfoDb::seqlevelsStyle(bt.genes) <- "UCSC"

    # plot 2
    output$overlap_plot <- renderPlot({
      #sel.genes <- subset(bt.genes, seqnames %in% paste0("chr", 1:1))
      sel.genes <- subset(bt.genes, seqnames %in% paste0("chr", input$chrII : input$chrII))
      # sel.genes <- subset(sel.genes, gene_biotype == "protein_coding")

      sel.genes <- subset(sel.genes, gene_biotype == input$biotype)
      #sel.cnvrs <- subset(cnvrs, seqnames %in% paste0("chr", 1:1))
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
      #sel.genes <- subset(bt.genes, seqnames %in% paste0("chr", 1:1))
      sel.genes <- subset(bt.genes, seqnames %in% paste0("chr", input$chrIII : input$chrIII))
      #sel.genes <- subset(sel.genes, gene_biotype == "protein_coding")

      sel.genes <- subset(sel.genes, gene_biotype == input$biotypeI)
      #sel.cnvrs <- subset(cnvrs, seqnames %in% paste0("chr", 1:1))
      sel.cnvrs <- subset(cnvrs, seqnames %in% paste0("chr", input$chrIII : input$chrIII))

      res <- suppressWarnings(overlapPermTest(A=sel.cnvrs, B=sel.genes, ntimes=100,
                                              genome="hg19", mask=NA, per.chromosome=TRUE, count.once=TRUE))
      summary(res[[1]]$permuted)
      #
      plot(res)
    })

  }

  # Run the application
  shinyApp(ui = ui, server = server)

}
