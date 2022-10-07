#load libraries
source("global.R", local = TRUE)
library(shiny)
library(shinydashboard)
library(DT)

# UI
ui <- fluidPage(
  
  dashboardPage(
    dashboardHeader(title = "CNVviz"),
    dashboardSidebar(
      sidebarMenu(
        menuItem("Import Data", tabName = "DataImport", icon = icon("dna")),
        menuItem("Recurrent Regions", tabName = "Recurrent_regions", icon = icon("dna")),
        menuItem("Overlap Analysis", tabName = "Overlap_analysis", icon = icon("dna")),
        menuItem("Overlap Permutation Test", tabName = "Overlap_permutation_test", icon = icon("dna"))
      )
    ),
    dashboardBody(use_waiter(),
                  tabItems(
                    
                    #### ------------------ Import data dashboard -------------------###############
                    tabItem(tabName = "DataImport", h2("Visualize CNVs data"),
                            #start fluidrow csv file load
                            fluidRow(
                              # start box right (contain input) ==> read csv file
                              box(status = "primary", width = 3,
                                  solidHeader = TRUE,
                                  collapsible = TRUE,
                                  title = "Import File",
                                  br(),
                                  "Import the variant data csv file here",
                                  br(),
                                  fileInput("input$csv_file", label = "Upload csv file"),
                                  br(),
                                  br(),
                                  br(),
                                  checkboxInput("checkbox", label = "Demo", value = FALSE),
                              ), # end box right
                              
                              # start box left (contain output) ==> ouput 
                              box(title = "Input file", width = 9, height = 600, 
                                  br(),
                                  DT::dataTableOutput("data_matrix_df")
                              ) # end box left 
                              
                            ) # end fluid row 1
                    ), #import tabitem

                    #### ------------------Recurrent regions -------------------###############
                    tabItem(tabName = "Recurrent_regions", h2("Identifying recurrent regions"),
                            box(title = "About", width = 3, height = 600, 
                                solidHeader = TRUE,
                                collapsible = TRUE,
                                "The The GISTIC method is used to identify those regions of the genome that are aberrant more often than would be expected by chance, with greater weight given to high amplitude events (high-level copy-number gains or homozygous deletions) that are less likely to represent random aberrations (Beroukhim et al., 2007).",
                                br(),
                                br(),
                                br(),
                                br(),
                                br(),
                                "Please define the significance threshold to filter recurrent CNVs.",
                                br(),
                                textInput("pvalue", label = h3("P value"), value = "0.05"),
                                br(),
                                "Please define the chromosome",
                                br(),
                                textInput("chr", label = h3("Chromosome"), value = "chr1")
                                ),
                            box(title = "Recurrent", width = 9, height = 800, 
                                br(),
                                "The plot summarize CNV regions, a valid UCSC genome assembly, and a chromosome of interest.",
                                br(),
                                plotOutput("recurrent_plot", height = 800))
                            ),
                    #### ------------------ Overlap analysis ------------------###############
                    tabItem(tabName = "Overlap_analysis", h2("Overlap analysis of CNVs with functional genomic regions"),
                            box(title = "About", width = 3, height = 200, 
                                solidHeader = TRUE,
                                collapsible = TRUE,
                                "It is of interest whether the resulting CNV regions overlap with functional genomic regions such as genes, promoters, or enhancers."),
                            box(title = "Overlap Plot", width = 9, height = 600, 
                                br(),
                                "Stacked barplots on the top and the right of the plot display the number of altered genes per sample and the number of altered samples per gene, respectively.",
                                br(),
                                plotOutput("overlap_plot", height = 600))
                    ),
                    #### ------------------ Overlap permutation test------------------###############
                    tabItem(tabName = "Overlap_permutation_test", h2("Overlap permutation test"),
                            #box(title = "About", width = 3, height = 200, 
                             #   solidHeader = TRUE,
                             #   collapsible = TRUE,
                             #   "..."),
                            box(title = "Overlap Plot", width = 9, height = 600, 
                                br(),
                                plotOutput("permutation_plot", height = 600))
                    )
                    
                  ) #tabitems
    ) #body
  ) #dashboard
  
) # fluidpage


# Server 
server <- function(input, output, session) {
  
  #    if(input$checkbox  == TRUE){
 # source("demo/cnv_viz.R", local = TRUE)
  #output$data_matrix_df <- DT::renderDataTable({
  #  req(calls)
  #  calls
  #}) 
#} else{
  #source("cnv_viz1.R", local = TRUE)
#}

  source("demo/cnv_viz.R", local = TRUE)

  ## read data csv file
  file_data <- reactive({
    file1 <- input$csv_file
    if(!is.null(file1)){read.csv(file1$datapath)}
  })
  
  output$data_matrix_df <- DT::renderDataTable({
    req(file_data())
    file_data()
  }) 
  
# using the input file
  source("cnv_viz1.R", local = TRUE)

}

# Run the application 
shinyApp(ui = ui, server = server)
