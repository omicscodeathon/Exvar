#load libraries
#source("global.R", local = TRUE)
library(shiny)
library(shinydashboard)
library(waiter)
library(DT)
library(echarts4r)
# UI
ui <- fluidPage(
  
  dashboardPage(
    dashboardHeader(title = "SNPviz"),
    dashboardSidebar(
      sidebarMenu(
        menuItem("Data Import", tabName = "DataImport", icon = icon("dna")),
        menuItem("SNP Data", tabName = "SNP", icon = icon("dna")),
        menuItem("Amino Acid Changes", tabName = "effect", icon = icon("dna"))
      )

    ),
    dashboardBody(use_waiter(),
                  tabItems(

                    #### ------------------ Import data dashboard -------------------###############
                    tabItem(tabName = "DataImport", 
                            #start fluidrow vcf file load
                            fluidRow(
                              # start box right (contain input) ==> read csv file
                              box(status = "primary", width = 3,
                                  solidHeader = TRUE,
                                  collapsible = TRUE,
                                  title = "Import File",
                                  br(),
                                  "Import the variant data vcf file here",
                                  br(),
                                  fileInput("vcf_file", label = "Upload vcf file"),
                              ), # end box right
                              
                              # start box left (contain output) ==> ouput 
                              box(title = "Visualize variantsfile", width = 9, height = 600, 
                                   br(),
                                  DT::dataTableOutput("data_matrix_vcf")
                              ) # end box left 
                              
                            ) # end fluid row 1
                    ), #import tabitem
                    
                    #### ------------------ SNP dashboard -------------------###############
                    tabItem(tabName = "SNP", 
                            fluidRow(width =12,
                                     box(width =6,
                                       "SNPs types Pie plot",  "pie of SNPs types",
                                       br(),
                                       echarts4rOutput("pie_snp_type", height = 500)
                                     ),
                                     box(width =6,
                                       "TiTv distribution",  "TiTv plot",
                                       br(),
                                       plotOutput("TiTv", height = 500)
                                     ),
                                     box(width =6,
                                       "Trinucleotide pattern plot",  "Trinucleotide pattern plot",
                                       br(),
                                       plotOutput("Tri_pattern", height = 500)
                                     ),
                                     box(width =6,
                                        "True and False variants count",
                                       br(),
                                       echarts4rOutput("Variantions", height = 500)
                                     ) 
                            ) #end fluid
                    ), #snp tabitem
                    
                    #### ------------------ AAC dashboard -------------------###############
                    tabItem(tabName = "effect",
                            fluidRow(width =12,
                                     box(title = "The Amino Acid Changes", width =12,
                                          br(),
                                          echarts4rOutput("Variantions", height = 500)
                                         )
                                     )
                    )
                  ) #tabitems
    ) #body
  ) #dashboard
  
) # fluidpage


# Server 
server <- function(input, output, session) {
  
  ## read vcf file
 file2_data <- reactive({
   file3 <- input$vcf_file
   if(!is.null(file3)){read.csv(file3$datapath)}
  })
  
  # start box 1 csv file  >>>>    how to skip the info lines !!!!!
  output$data_matrix_vcf <- DT::renderDataTable({
    req(file2_data())
    file2_data()
  }) 
  ####---------------------- filter   ---------------------------####
  #source("server/variant_filter.R", local = TRUE)
  ####---------------------- snp types visualization   ---------------------------####
  source("server/snp_viz.R", local = TRUE)
  # 
  output$data_matrix_vcf <- DT::renderDataTable({
    req(snp_data)
    snp_data
  }) 
  ####---------------------- snp effect   ---------------------------####
  source("server/snp_effect.R", local = TRUE)
  ####---------------------- Amino Acids changes  ---------------------------####
  source("server/aach.R", local = TRUE)
}

# Run the application 
shinyApp(ui = ui, server = server)
