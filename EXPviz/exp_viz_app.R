#load libraries
#source("global.R", local = TRUE)

# UI
ui <- fluidPage(
  
  dashboardPage(
    dashboardHeader(title = "EXPviz"),
    dashboardSidebar(
      sidebarMenu(
        menuItem("Import Data", tabName = "ImportData", icon = icon("dna")),
        menuItem("Count Data", tabName = "CountData", icon = icon("dna")),
        menuItem("Expression", tabName = "Expression", icon = icon("dna")),
        menuItem("Ontology", tabName = "Ontology", icon = icon("dna"))
      )

    ),
    dashboardBody(use_waiter(),
                  tabItems(
                    
                    #### ------------------ Import data dashboard -------------------###############
                    tabItem(tabName = "ImportData",
                            #start fluidrow vcf file load
                            fluidRow(
                              # start box right (contain input) ==> read csv file
                              box(status = "primary", width = 2,
                                  solidHeader = TRUE,
                                  collapsible = TRUE,
                                  title = "Import count File",
                                  br(),
                                  "Import the count csv file here",
                                  br(),
                                  fileInput("csv_file", label = "count.csv")
                              ), # end box right
                              box(status = "primary", width = 2,
                                  solidHeader = TRUE,
                                  collapsible = TRUE,
                                  title = "Import metadat File",
                                  br(),
                                  "Import the metadata txt file here",
                                  br(),
                                  fileInput("metadata_file", label = "metadata.txt")
                              ), # end box right                              
                              # start box left (contain output)
                              box(title = "Visualize csv file", width = 8, height = 700, 
                                   br(),
                                  DT::dataTableOutput("data_matrix_csv")
                              ) # end box left 
                              
                            ) # end fluid row 1
                    ), #tab item tabitem
                    
                    #### ------------------ count -------------------###############
                    tabItem(tabName = "CountData",
                    fluidRow(
                      # start bow right
                      tabBox(title = "Read Count before log2 transformation", width = 6,
                             ### multiple plots in same box
                             selected = "Read_Count1",
                             
                             tabPanel("Read_Count1", "Before log2 transformation", 
                                      br(),
                                      plotOutput("read_count1", height = 300),
                                      plotOutput("read_count2", height = 300))),
                      tabBox(title = "Read Count after log2 transformation", width = 6,
                             ### multiple plots in same box
                             selected = "Read_Count2",
                             
                             tabPanel("Read_Count2", "After log2 transformation",
                                      br(),
                                      plotOutput("log_transformation1", height = 300),
                                      plotOutput("log_transformation2", height = 300)),
                             tabPanel("density plot ",  "density plot ",
                                      br(),
                                      plotOutput("density_plot")),
                             
                             tabPanel("scatter plot ",  "scatter plot ",
                                      br(),
                                      textInput("sam1", label = h3("sample ID"), value = "eg: SRR11426825"),
                                      br(),
                                      textInput("sam1", label = h3("sample ID"), value = "eg: SRR11426823"),
                                      br(),
                                      plotOutput("scatter_plot"))
                    ))),
                    #### ------------------ expression -------------------###############
                    tabItem(tabName = "Expression",
                            
                            fluidRow(

                              box(status = "info", width = 6,
                                  solidHeader = TRUE,
                                  collapsible = TRUE,
                                  title = "expression",
                                  br(),
                                  textInput("pval", label = h3("P value"), value = "0.05"),
                                  br(),
                                  textInput("log", label = h3("LogFC"), value = "2")),
                              br(),   
                               box(
                                 "PCA",  "PCA plot", width = 6,
                                  br(),
                                  plotOutput("pca_plot", height = 500)
                               ),
                              br(),
                               box(
                                 "volcano", "Volcano plot",width = 6,
                                 br(),
                                 plotOutput("volcano_plot", height = 500)
                               ),
                              br(),
                               box(
                                 "MA", "MA plot",width = 6,
                                 br(),
                                 plotOutput("ma_plot", height = 500)
                               )
                               # end box 
                            ),
                            
                            fluidRow(width =12,
                                     tabBox(title = "Heatmaps", width =12,
                                            
                                            selected = "Heatmap",
                                         
                                            tabPanel("Heatmap", 
                                                     br(),
                                                     "Heatmap of all expressed genes",
                                                     br(),
                                                     "Note :The heatmap need 2 minutes to load.",
                                                     plotOutput("heatmap", height = 500)),
                                            tabPanel("HeatmapDEGs",
                                                     br(),
                                                     "Heatmap of differentially expressed genes only",
                                                     br(),
                                                     "Note :The heatmap need 2 minutes to load.",
                                                     plotOutput("heatmap_degs", height = 500))
                                            
                                     ) # end tabbox
                            ) # 
                            ),
                    
                    
                    #### ------------------ ontology -------------------###############
                    tabItem(tabName = "Ontology",
                            fluidRow(width =12,
                            box(status = "info", width = 3,
                                solidHeader = TRUE,
                                collapsible = TRUE,
                                title = "Filters",
                                br(),
                                #selectbox
                                "Genes to analyze",
                                "Note: DEGs= Differentially expressed genes",
                                selectInput("select", label = h3("Choose gene group"), 
                                            choices = list("all DEGs" = 1, "Upregulated genes" = 2, "Downregulated genes" = 3), 
                                            selected = 1),
                                br(),
                                br(),
                                br(),
                                br(),
                                br(),
                                textInput("show", label = h3("Top Ontologies"), value = "20")),
                            tabBox(title = "Ontologies", width =9,
                                   
                                   fluidRow(width =12,
                                            #BP
                                            tabBox(title = "Biological Processes", width =12,
                                                   ### multiple plots in same box
                                                   selected = "Braplot1",
                                                   tabPanel("Barplot1",  "BP Barplot",
                                                            br(),
                                                            plotOutput("BP_barplot", height = 500)),
                                                   
                                                   tabPanel("Dotplot", "BP dotplot",
                                                            br(),
                                                            plotOutput("BP_dotplot", height = 500))
                                            ), # end tabbox
                                            #MF
                                            tabBox(title = "Molecular Functions", width =12,
                                                   ### multiple plots in same box
                                                   selected = "Braplot2",
                                                   tabPanel("Barplot2",  "MF Barplot",
                                                            br(),
                                                            plotOutput("MF_barplot", height = 500)),
                                                   
                                                   tabPanel("Dotplot", "MF dotplot",
                                                            br(),
                                                            plotOutput("MF_dotplot", height = 500))
                                            ), # end tabbox
                                            #CC
                                            tabBox(title = "Cellular Compenent", width =12,
                                                   ### multiple plots in same box
                                                   selected = "Braplot3",
                                                   tabPanel("Barplot3",  "CC Barplot",
                                                            br(),
                                                            plotOutput("Cc_barplot", height = 500)),
                                                   
                                                   tabPanel("Dotplot", "Cc dotplot",
                                                            br(),
                                                            plotOutput("CC_dotplot", height = 500))
                                            ) # end tabbox
                                   ) #end see vcf fluid
                            )
                            
                            )#fuidrow1
                            )#tabitem

                    
                  ) #tabitems
    ) #body
  ) #dashboard
  
) # fluidpage


# Server 
server <- function(input, output, session) {
  
  #### ------------------ expression dashboard -------------------###############
  ## read count data csv file
  file_data <- reactive({
    file1 <- input$csv_file
    if(!is.null(file1)){read.csv(file1$datapath)}
  })
  
  # start box 1 csv file
  output$data_matrix_csv <- DT::renderDataTable({
    req(file_data())
    file_data()
  }) 
  ## read metadat data file
  meta_file <- reactive({
    file2 <- input$metadata_file
    if(!is.null(file2)){read.csv(file2$datapath)}
  })
  
  # expression interface
  ####---------------------- viz count data   ---------------------------####
 source("server/count_viz.R", local = TRUE)
  ####---------------------- viz exp data   ---------------------------####
 source("server/exp_viz.R", local = TRUE)
  ####---------------------- Ontology analysis   ---------------------------####
  source("server/ontology.R", local = TRUE)
  ####---------------------- Ontology analysis   ---------------------------####
  #source("server/pathway.R", local = TRUE)
}

# Run the application 
shinyApp(ui = ui, server = server)
