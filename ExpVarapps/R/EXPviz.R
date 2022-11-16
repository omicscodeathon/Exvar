#' EXPviz
#'
#' @param countdata path to the count data csv file
#' @param metadata path to the metadata csv file
#'
#' @return
#' @export
#'
#' @examples
EXPviz <- function(countdata,metadata) {

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
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(enrichplot)
  library(shinydashboard)

  countdata <- as_data_frame(read_csv(countdata))
  metadata <- as.data.frame(read_table(metadata))
  ui <- bootstrapPage(
    navbarPage(theme = shinytheme("flatly"), collapsible = TRUE,
               HTML('<a style="text-decoration:none;cursor:default;color:#FFFFFF;" class="active" href="#">EXPviz</a>'), id="nav",
               windowTitle = "EXPviz",
               tabPanel("Count Data",
                        tabBox(title = "Before log2 transformation", width = 6,
                               selected = "Read Count",
                               tabPanel("Read Count",
                                        plotOutput("read_count1", height = 300)),
                               tabPanel("Histogram",
                                        plotOutput("read_count2", height = 300))),
                        tabBox(title = "After log2 transformation", width = 6,
                               selected = "Read Count",

                               tabPanel("Read count",
                                        plotOutput("log_transformation1", height = 300)),
                               tabPanel("Read count",
                                        plotOutput("log_transformation2", height = 300)),
                               tabPanel("density plot ",
                                        plotOutput("density_plot"))
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
                            span(shiny::tags$i(h6("Define the LogFC to identify the differentially expressed genes")), style="color:#045a8d"),
                            br(),
                            sliderInput("log",
                                        "Log value:",
                                        min = 0,
                                        max = 5,
                                        value=2
                            ),
                            br(),
                            varSelectInput("factor", "DESeq Factor:", metadata)
                          ),
                          mainPanel(
                            tabsetPanel(
                              tabPanel("Volcano Plot", downloadButton("downloadvplot", "Download"),
                                       plotOutput("volcano_plot", height = 500)),
                              tabPanel("MA Plot", downloadButton("downloadmplot", "Download"),
                                       plotOutput("ma_plot", height = 500)),
                              tabPanel("PCA Plot", downloadButton("downloadpplot", "Download"),
                                       plotOutput("pca_plot", height = 500)),
                              tabPanel("Heatmap", downloadButton("downloadhplot", "Download"),
                                       plotOutput("heatmap", height = 500))
                            )
                          )
                        )),
               tabPanel("Ontology",

                        sidebarLayout(
                          sidebarPanel(

                            span(shiny::tags$i(h6("Define the gene group")), style="color:#045a8d"),
                            br(),
                            pickerInput("select", "Gene Group:",
                                        choices = c("All Differentially Expresssed Genes", "Upregulated genes", "Downregulated genes" ),
                                        selected = c("All Differentially Expresssed Genes"),
                                        multiple = FALSE),
                            br(),
                            span(shiny::tags$i(h6("Top Ontologies")), style="color:#045a8d"),
                            br(),
                            textInput("show", label="Top ontologies",value = "20")
                          ),
                          mainPanel(
                            tabsetPanel(
                              tabPanel("Biological Processes", downloadButton("downloadBPplot", "Download"),
                                       plotOutput("BP_dotplot", height = 500)),
                              tabPanel("Cellular Compenant", downloadButton("downloadCCplot", "Download"),
                                       plotOutput("CC_dotplot", height = 500)),
                              tabPanel("Molecular Function", downloadButton("downloadMFplot", "Download"),
                                       plotOutput("MF_dotplot", height = 500))
                            )
                          )
                        ))
    ))

  server <- function(input, output) {

    countdata = countData
    rownames(countdata) = countData$gene   #first column name
    countData = countData[,-c(1)]


    output$read_count1 <- renderPlot({
      par(mar=c(8,4,4,1)+0.1)
      barplot( colSums(countData)/1e6, las=3, main="Total read counts (millions)", ylab="Total read counts in millions")

    })
    ## error
    # output$read_count2 <- renderPlot({
    #   hist(countData[,1], br=200, xlab="Number of Reads Counts per Feature", main="Histogram of Read Counts")

    # })
    # log transformation + log transformed read count plot
    logCountData = log2(1+countData)
    par(mfrow = c(1, 2), mar=c(8,4,4,1))  # two columns

    #
    output$log_transformation1 <- renderPlot({
      hist(logCountData[,1], main="Histogram of Log Read Counts", xlab="Log transformed counts")
    })
    output$log_transformation2 <- renderPlot({
      boxplot(logCountData,las=3, main="Boxplot of Log Read Counts")
    })
    # density plot
    output$density_plot <- renderPlot({
      logCountData = log2(1+countData)
      par(mfrow = c(1, 2), mar=c(8,4,4,1))  # two columns   change margins
      x <- logCountData
      myColors = rainbow(dim(x)[2])
      plot(density(x[,1]),col = myColors[1], lwd=2,
           xlab="Expresson values", ylab="Density", main= "Distribution of transformed data",
           ylim=c(0, max(density(x[,1])$y)+.02 ) )
      for( i in 2:dim(x)[2] )
        lines(density(x[,i]),col=myColors[i], lwd=2)
      legend("topright", cex=1.1,colnames(x), lty=rep(1,dim(x)[2]), col=myColors )
    })
    #### expression

    coldata <- metadata
    coldata %>% remove_rownames %>% column_to_rownames(var="sample")
    countData = as.data.frame(countdata)
    rownames(countData) = countData$gene   #first column name
    countData = countData[,-c(1)]
    counts_data <- countData
    # creat deseq dataset object  ## add button Run DE
    reactive_df <- reactive({
      dds <- DESeqDataSetFromMatrix(countData = counts_data,
                                    colData = coldata,
                                    design = ~ input$factor)

      dds
    })
    # output$data <- renderTable({
    # factor = ~ factor1  #input
    #  dds <- DESeqDataSetFromMatrix(countData = counts_data,
    #                                colData = coldata,
    #                                design = ~ input$factor)
    #
    #   dds
    # })

    # pre-filtering: removing rows with low gene counts
    # keeping rows that have at least 10 reads total
    keep <- rowSums(counts(dds)) >= 10  # reads number could be customized
    dds <- dds[keep,]

    dds
    # Step 3: Run DESeq ----------------------
    dds <- DESeq(dds)
    expression_data  <- results(dds)

    # MA plot
    output$ma_plot <- renderPlot({
      DESeq2::plotMA(expression_data,  ylim = c(-5, 5))
    })

    # identify DEGss

    # volcano plot
    output$volcano_plot <- renderPlot({

      expression = as.data.frame(expression_data)
      expression$diffexpressed <- "Not Significant"
      expression$diffexpressed[expression$log2FoldChange  > input$log  & expression$pvalue < input$pval ] <- "Upregulated genes"
      expression$diffexpressed[expression$log2FoldChange  < -(input$log)  & expression$pvalue < input$pval ] <- "Downregulated genes"

      ggplot(data=expression, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed)) + geom_point() + theme_minimal()

    })

    #pca
    output$pca_plot <- renderPlot({
      rld = rlog(dds)
      # factorpca <- input$factor
      plotPCA(rld, intgroup= input$factor )

    })
    # heatmap all genes

    output$heatmap <- renderPlot({
      mat <-counts(dds, normalized=T)[rownames(expression),]
      mat.z <- t(apply(mat,1,scale))
      cold = coldata
      rownames(cold) <- cold[,1]
      colnames(mat.z)<-rownames(cold)
      mat.z

      Heatmap(mat.z, cluster_rows = T, cluster_columns = T, column_labels = colnames(mat.z), name="Z-score",
      )

    })

    ########
    #filter upregulated genes
    expression <- tibble::rownames_to_column(expression, "ENSEMBL")
    upgene <- as.data.frame(expression[expression$diffexpressed == "Upregulated genes",])
    upgenes_list <- upgene$ENSEMBL

    downgene <- as.data.frame(expression[expression$diffexpressed == "Downregulated genes",])
    downgenes_list <- downgene$ENSEMBL

    all_degs <- as.data.frame(expression[expression$diffexpressed != "Not Significant",])
    all_degs_list <- downgene$ENSEMBL
    #

    output$BP_dotplot<- renderPlot({
      if (input$select == "All Differentially Expresssed Genes"){
        group_gene <- all_degs_list
      } else if (input$select == "Upregulated genes"){
        group_gene <- upgenes_list
      } else {
        group_gene <- downgenes_list
      }


      BP_up_go <- enrichGO(gene = group_gene, OrgDb="org.Hs.eg.db", keyType="ENSEMBL", ont="BP")
      BP_up_table <- as.data.frame(BP_up_go)

      dotplot(BP_up_go, showCategory=30)
    })



    output$MF_dotplot<- renderPlot({
      if (input$select == "All Differentially Expresssed Genes"){
        group_gene <- all_degs_list
      } else if (input$select == "Upregulated genes"){
        group_gene <- upgenes_list
      } else {
        group_gene <- downgenes_list
      }

      MF_up_go <- enrichGO(gene = group_gene, OrgDb="org.Hs.eg.db", keyType="ENSEMBL", ont="MF")
      dotplot(BP_up_go, showCategory=30)
    })


    output$CC_dotplot<- renderPlot({
      if (input$select == "All Differentially Expresssed Genes"){
        group_gene <- all_degs_list
      } else if (input$select == "Upregulated genes"){
        group_gene <- upgenes_list
      } else {
        group_gene <- downgenes_list
      }

      CC_up_go <- enrichGO(gene = group_gene, OrgDb="org.Hs.eg.db", keyType="ENSEMBL", ont="CC")
      dotplot(CC_up_go, showCategory=30)
    })


  }

  # Run the application
  shinyApp(ui = ui, server = server)

}
