#'  vizexp
#'
#' @param genecount Path to the gene count file
#' @param metadata Path to the  metadata file
#' @param species Species name ( "Homo Sapiens", "Mus Musculus","Arabidopsis Thaliana" ,
#' "Drosophila Melanogaster", "Danio rerio", "Rattus norvegicus", or "Saccharomyces Cerevisiae" )
#'
#'
vizexp <- function( genecount, metadata, species) {

 #load requirements
  library(shinydashboard)
  library(shiny)
  library(shinythemes)
  library(shinyWidgets)
  library(readr)
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(clusterProfiler)
  library(AnnotationDbi)
  library(enrichplot)
  library(ggnewscale)
  library(tibble)
  library(readxl)
  library(shinycssloaders)

  #specie
  if (species == "Homo Sapiens"){
    library(org.Hs.eg.db)
    specie1 <- "org.Hs.eg.db"
  } else if (species == "Mus Musculus"){
    library(org.Mm.eg.db)
    specie1 <-"org.Mm.eg.db"
  } else if (species == "Arabidopsis Thaliana"){
    library(org.At.tair.db)
    specie1 <-"org.At.tair.db"
  } else if (species == "Drosophila Melanogaster"){
    library(org.Dm.eg.db)
    specie1 <-"org.Dm.eg.db"
  } else if (species == "Danio rerio"){
    library(org.Dr.eg.db)
    specie1 <-"org.Dr.eg.db"
  } else if (species == "Rattus_norvegicus"){
    library(org.Rn.eg.db)
    specie1 <-"org.Rn.eg.db"
  } else if (species == "Saccharomyces Cerevisiae"){
    library(org.Sc.sgd.db)
    specie1 <-"org.Sc.sgd.db"
  }

  #prepare files
  #count data
  countdata <- read_csv(genecount)
  countdata = as.data.frame(countdata)
  # Remove duplicates based on gene column
  countdata <- countdata[!duplicated(countdata$gene), ]
  countdata<- na.omit(countdata)

  rownames(countdata) = countdata$gene   #first column name
  countdata  = countdata [,-c(1)]

  #metadata
  metadata <- read_excel("metadata.xlsx")
  coldata <-as.data.frame( na.omit(metadata))
  rownames(coldata) = coldata$Sample

  # #run app
  ui <- bootstrapPage(
    navbarPage(theme = shinytheme("flatly"), collapsible = TRUE,
               HTML('<a style="text-decoration:none;cursor:default;color:#FFFFFF;" class="active" href="#">EXPviz</a>'), id="nav",
               windowTitle = "EXPviz",
               tabPanel("Count Data",
                        tabBox(title = "Before log2 transformation", width = 6,
                               box(shinycssloaders::withSpinner(plotOutput("read_count1")))),

                        tabBox(title = "After log2 transformation", width = 6,
                               tabPanel("Read count",
                                        shinycssloaders::withSpinner(plotOutput("log_transformation2"))),
                               tabPanel("density plot ",
                                        shinycssloaders::withSpinner(plotOutput("density_plot")))
                        )),
               tabPanel("Expression Analysis",
                        sidebarLayout(
                          sidebarPanel(

                            span(shiny::tags$i(h6("Define the p-value to identify the differentially expressed genes")), style="color:#045a8d"),
                            br(),
                            sliderInput("pval",
                                        "P_value:",
                                        min = 0,
                                        max = 1,
                                        value=0.05
                            ),
                            br(),
                            span(shiny::tags$i(h6("Define the LogFC value to identify the differentially expressed genes")), style="color:#045a8d"),
                            br(),
                            sliderInput("log",
                                        "Log value:",
                                        min = 0,
                                        max = 5,
                                        value=2
                            )
                          ),
                          mainPanel(
                            tabsetPanel(
                              tabPanel("MA Plot",
                                       shinycssloaders::withSpinner(plotOutput("ma_plot", height = 500))

                              ),
                              tabPanel("PCA Plot",
                                       shinycssloaders::withSpinner(plotOutput("pca_plot", height = 500))

                              ),
                              tabPanel("Volcano Plot",
                                       shinycssloaders::withSpinner(plotOutput("volcano_plot", height = 500))
                                       )
                            )
                          )
                        )),
               tabPanel("Ontology",

                        sidebarLayout(
                          sidebarPanel(
                            pickerInput("select", label= span(shiny::tags$i(h6("Define the gene group")), style="color:#045a8d"),

                                        choices = c("All Differentially Expresssed Genes", "Upregulated genes", "Downregulated genes" ),
                                        selected = c("All Differentially Expresssed Genes"),
                                        multiple = FALSE),
                            br(),
                            textInput("show", label=  span(shiny::tags$i(h6("Top Ontologies")), style="color:#045a8d")
                                      ,value = "20"),
                            br(),
                            pickerInput("plottype", label= span(shiny::tags$i(h6("Define the Plot type")), style="color:#045a8d"),

                                        choices = c("Dotplot","Barplot",  "Cnetplot" ),
                                        multiple = FALSE),
                            br(),
                            pickerInput("GOtype", label= span(shiny::tags$i(h6("Define the Ontology category")), style="color:#045a8d"),

                                        choices = c("Biological Processes","Cellular Compenant",  "Molecular Function" ),
                                        multiple = FALSE)
                          ),
                          mainPanel(
                            shinycssloaders::withSpinner(plotOutput("ontology_plot"))

                          )
                        )
               )
    ))

  server <- function(input, output) {

    output$read_count1 <- renderPlot({
      par(mar=c(8,4,4,1)+0.1)
      barplot( colSums(countdata)/1e6, las=3, main="Total read counts (millions)", ylab="Total read counts in millions")
    })

    # log transformation
    logData = log2(1+countdata)
    par(mfrow = c(1, 2), mar=c(8,4,4,1))
    #
    output$log_transformation2 <- renderPlot({
      boxplot(logData,las=3, main="Boxplot of Log Read Counts")
    })
    # density plot
    output$density_plot <- renderPlot({
      myColors = rainbow(dim(logData)[2])
      plot(density(logData[,1]),col = myColors[1], lwd=2,
           xlab="Expresson values", ylab="Density", main= "Distribution of transformed data",
           ylim=c(0, max(density(logData[,1])$y)+.02 ) )
      for( i in 2:dim(logData)[2] )
        lines(density(logData[,i]),col=myColors[i], lwd=2)
      legend("topright", cex=1.1,colnames(logData), lty=rep(1,dim(logData)[2]), col=myColors )
    })

    #expression analysis
    dds <- DESeqDataSetFromMatrix(countData = countdata,
                                  colData = coldata,
                                  design = ~Factor)
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]
    dds <- DESeq(dds)
    expression_data  <-  results(dds)
    expression <- as.data.frame(expression_data)
    expression <- rownames_to_column(expression, "gene")
    expression$diffexpressed <- "Not Significant"
    # MA plot
    output$ma_plot <- renderPlot({
      DESeq2::plotMA(expression_data,  ylim = c(-5, 5))
    })

    #pca plot
    output$pca_plot <- renderPlot({
      rld = rlog(dds)
      plotPCA(rld, intgroup= "Factor")
    })

    # volcano plot
    output$volcano_plot <- renderPlot({
      #define degs
      expression$diffexpressed[expression$log2FoldChange  > input$log  & expression$pvalue < input$pval ] <- "Upregulated genes"
      expression$diffexpressed[expression$log2FoldChange  < -(input$log)  & expression$pvalue < input$pval ] <- "Downregulated genes"
      #plot
      ggplot(data=expression, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed)) + geom_point() + theme_minimal()

    })

    #ontolgy analysis
    output$ontology_plot<- renderPlot({
      #define degs
      expression$diffexpressed[expression$log2FoldChange  > input$log  & expression$pvalue < input$pval ] <- "Upregulated genes"
      expression$diffexpressed[expression$log2FoldChange  < -(input$log)  & expression$pvalue < input$pval ] <- "Downregulated genes"
      #gene group
      if (input$select == "All Differentially Expresssed Genes"){
        all_degs <- as.data.frame(expression[expression$diffexpressed != "Not Significant",])
        group_gene <- all_degs$gene
      } else if (input$select == "Upregulated genes"){
        upgene <- as.data.frame(expression[expression$diffexpressed == "Upregulated genes",])
        group_gene <- upgene$gene
      } else {
        downgene <- as.data.frame(expression[expression$diffexpressed == "Downregulated genes",])
        group_gene <- downgene$gene
      }
      #ontology type

      if (input$GOtype == "Biological Processes"){
        GO_type <- "BP"
      }else if (input$GOtype == "Cellular Compenant"){
        GO_type <- "CC"
      }else{
        GO_type <- "MF"
      }
      # analysis
      ontology <- enrichGO(gene = group_gene, OrgDb=specie1, keyType="SYMBOL", ont=GO_type)
      s_ontology <- simplify(ontology)
      #plot
      if (input$plottype == "Dotplot"){
        dotplot(ontology, showCategory=as.numeric(input$show))
      }else if (input$plottype == "Barplot"){
        barplot(ontology, showCategory=as.numeric(input$show))
      }else{
        cnetplot(s_ontology, foldChange=group_gene, circular = TRUE, colorEdge = TRUE)
      }

    })

  }
  # Run the application
  shinyApp(ui = ui, server = server)

}
