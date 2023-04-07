#' vizexp
#' This function is a shiny app to visualize gene expression data.
#' @param genecount Path to the gene count file
#' @param metadata Path to the  metadata file
#' @export
#'
vizexp <- function( genecount, metadata) {
  #request specie
  cat(paste0("These are the species currently supported by Exvar: \n",
             "[1] Homo sapiens \n",
             "[2] Mus musculus \n",
             "[3] Arabidopsis thaliana \n",
             "[4] Drosophila melanogaster \n",
             "[5] Danio rerio \n",
             "[6] Rattus norvegicus \n",
             "[7] Saccharomyces cerevisiae \n",
             "[8] Caenorhabditis elegans \n"))
  species <- readline("These are the species currently supported by Exvar, type the number corresponding to the target specie: ")
  #  load requirements
  library(shinydashboard) 
  library(shiny)
  library(shinythemes)
  library(shinyWidgets)
  library(shinycssloaders)
  library(readr)
  library(DESeq2)
  library(dplyr) 
  library(ggplot2) 
  library(tibble)
  library(readxl)
  library(clusterProfiler)
  library(enrichplot)
  # species specific requirements
  switch(species,
         "1"={
           library(org.Hs.eg.db)
           specie1 <- "org.Hs.eg.db"
         },
         "2"={
           library(org.Mm.eg.db)
           specie1 <-"org.Mm.eg.db"
         },
         "3"={
           library(org.At.tair.db)
           specie1 <-"org.At.tair.db"
         },
         "4"={
           library(org.Dm.eg.db)
           specie1 <-"org.Dm.eg.db"
         },
         "5"={
           library(org.Dr.eg.db)
           specie1 <-"org.Dr.eg.db"
         },
         "6"={
           library(org.Rn.eg.db)
           specie1 <-"org.Rn.eg.db"
         },
         "7"={
           library(org.Sc.sgd.db)
           specie1 <-"org.Sc.sgd.db"
         },
         "8"={
           library(org.Ce.eg.db)
           specie1 <-"org.Ce.eg.db"
         }
  )
  ## read input files + process files      
  # Note the genecount file must contain only one column of the genes symbol named "genes"
  # Remove duplicates based on gene column
  countdata <- read.csv(genecount) %>%
    distinct(gene, .keep_all = TRUE) %>%
    na.omit() %>%
    as_tibble() %>%
    column_to_rownames("gene") %>%
    as.data.frame() %>%
    select(-1)
  
  #metadata = coldata
  # Note the metadata file must contain a column contained the samples ids named "sample"
  # and a column of the samples' phenotype named "factor"
  coldata <- read_excel(metadata) %>%
    na.omit() %>%
    column_to_rownames("sample")
  ##run app
  ui <-   navbarPage(
    theme = shinytheme("flatly"),
    HTML('<a style="text-decoration:none;
               cursor:default;
                    color:white;
                    " class="active" href="#">Visualize Gene Expression data</a>'),
    windowTitle ="expviz",
    tabPanel("Count Data",
             span(shiny::tags$i(h2("Visualize gene count data : ")), style="color:#045a8d"),
             tabBox(title = "Before log2 transformation", width = 6,
                    box(shinycssloaders::withSpinner(plotOutput("read_count1")))),
             
             tabBox(title = "After log2 transformation", width = 6,
                    tabPanel("Read count",
                             shinycssloaders::withSpinner(plotOutput("log_transformation2"))),
                    tabPanel("density plot ",
                             shinycssloaders::withSpinner(plotOutput("density_plot")))
             )),
    tabPanel("Expression Analysis",
             span(shiny::tags$i(h2("Perform gene expression analysis :")), style="color:#045a8d"),
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
    tabPanel("Ontology Analysis",
             span(shiny::tags$i(h2("Perform gene ontology analysis")), style="color:#045a8d"),
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
  )
  ## server
  server <- function(input, output,session) {
    
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
                                  design = ~factor)
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
      plotPCA(rld, intgroup= "factor")
    })
    
    # volcano plot
    output$volcano_plot <- renderPlot({
      #define degs
      expression$diffexpressed[expression$log2FoldChange  > input$log  & expression$pvalue < input$pval ] <- "Upregulated genes"
      expression$diffexpressed[expression$log2FoldChange  < -(input$log)  & expression$pvalue < input$pval ] <- "Downregulated genes"
      #plot
      ggplot(data=expression, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed)) +
        geom_point() + theme_minimal() + theme(text = element_text(size = 15))  
      
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
  
  # run app
  shinyApp(ui = ui, server = server)
}
