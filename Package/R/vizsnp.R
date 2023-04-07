#'  vizsnp
#'
#' @param dir Path contain two folders named "Control" and "patient" containing the VCF files.
#' @export
#'
vizsnp <- function( dir) {
  #load requirement
  library(shiny)
  library(shinydashboard)
  library(shinythemes)
  library(ggplot2)
  library(data.table)
  library(fs)
  library(dplyr)
  library(purrr)
  library(DT)
  
  ###function input : dir
  setwd(dir)
  #read files
  file_paths_control <- fs::dir_ls("control")
  file_paths_patient <- fs::dir_ls("patient")
  
  ###prepare data
  #part I : controls
  #loop
  file_contenets_control <- list()
  file_contenets_control <-file_paths_control %>%
    map(function(path){
      fread(path)
    })
  #add each list item into the environment variable as a tibble.
  # and set each list item to an environment variable
  file_contenets_control%>% map(as_tibble) %>%
    list2env(envir = .GlobalEnv)
  # To join all files together into one data frame7
  comined_control <- file_contenets_control %>% map(as_tibble) %>%
    reduce(full_join)
  # delete NAs
  # result : keep commun between all samples
  comined_control <- na.omit(comined_control)
  ## identify SNPs
  # Transition (Ti)
  ti <- c("A>G","G>A","C>T","T>C")
  # Transveersion (Tv)
  tv <- c("A>T","A>C","G>T","G>C","C>A","C>G","T>A","T>G")
  #
  comined_control <- comined_control %>%
    mutate(nuSub = paste0(REF, ">", ALT),
           TiTv = case_when(
             nuSub %in% ti ~ "Ti",
             nuSub %in% tv ~ "Tv",
             TRUE ~ NA_character_
           )) %>%
    rename(snp_controls = nuSub, TiTv_controls = TiTv)

  ### Part II : patients
  #path
  file_contenets_patient <- list()
  file_contenets_patient <-file_paths_patient %>%
    map(function(path){
      fread(path)
    })
  #add each list item into the environment variable as a tibble.
  # and set each list item to an environment variable
  file_contenets_patient%>% map(as_tibble) %>%
    list2env(envir = .GlobalEnv)
  # To join all files together into one data frame7
  comined_patient <-file_contenets_patient %>% map(as_tibble) %>%
    reduce(full_join)
  # delete NAs
  # result : keep commun between all samples
  comined_patient <- na.omit(comined_patient)
  ## identify snps
  comined_patient <- comined_patient %>%
    mutate(nuSub = paste0(REF, ">", ALT),
           TiTv = case_when(
             nuSub %in% ti ~ "Ti",
             nuSub %in% tv ~ "Tv",
             TRUE ~ NA_character_
           )) %>%
    rename(snp_patient = nuSub, TiTv_patient = TiTv)

  # Part III: compare groups

  compare_group <- merge(comined_patient, comined_control, by = c("#CHROM", "POS"), all = TRUE) %>%
    mutate(across(everything(), ~ifelse(is.na(.), "no", .))) %>%
    mutate(compare = case_when(
      snp_controls != "no" & snp_patient == "no" ~ "deletion",
      snp_controls == "no" & snp_patient != "no" ~ "addition",
      snp_controls == snp_patient ~ "same",
      snp_controls != snp_patient & snp_controls != "no" & snp_patient != "no" ~ "different",
      TRUE ~ NA_character_
    ))
  
  # what to keep ?
  uncommun <- compare_group[compare_group$compare == "different",]
  commun <- compare_group[compare_group$compare == "same",]

  #
  addition <- compare_group %>%
    filter(compare == "addition") %>%
    select("#CHROM", POS, REF.x, ALT.x, TiTv_patient) %>%
    rename(Chromosome = "#CHROM", `SNP position` = POS, REF = REF.x, ALT = ALT.x, TiTv = TiTv_patient)

  #
  deletion <- compare_group %>%
    filter(compare == "deletion") %>%
    select("#CHROM", POS, REF.y, ALT.y, TiTv_controls) %>%
    rename(Chromosome = "#CHROM", `SNP position` = POS, REF = REF.y, ALT = ALT.y, TiTv = TiTv_controls)
  
  # final dataframe 
  uncommun_simple <- uncommun %>%
    select(`#CHROM`, POS, ID.x, REF.x, ALT.x, ALT.y) %>%
    rename(CHROM = `#CHROM`, ID = ID.x, REF = REF.x, ALT_patient = ALT.x, ALT_control = ALT.y)
  #app
  ui <- fluidPage(
    navbarPage(theme = shinytheme("flatly"), collapsible = TRUE,
               HTML('<a style="text-decoration:none;cursor:default;color:#FFFFFF;" class="active" href="#">SNPviz</a>'),
               id="nav",
               windowTitle = "VizSNP",
               tabPanel("Single-nucleotide polymorphism",
                        h2("Compare SNPs between two groups"),
                        box( width = 6,
                             downloadButton("downloadData", "Download csv file"),
                             dataTableOutput('table_compare')
                        ),
                        tabBox(width = 6,
                               tabPanel("SNP distribution",
                                        shinycssloaders::withSpinner( plotOutput("snp_patient"))
                               ),
                               tabPanel("TiTv distribution",
                                        shinycssloaders::withSpinner( plotOutput("TiTv_patient"))
                               )
                        )),#tabpane
               tabPanel("Insertion/Deletion Polymorphism",
                        
                        box( width = 6,
                             h2( "These SNPs exist in the patient but not in the control group"),
                             br(),
                             br(),
                             downloadButton("downloadins", "Download csv file"),
                             br(),
                             br(),
                             dataTableOutput('insertion')
                        ),
                        box( width = 6,
                             h2( "These SNPs exist in the control but not in the patient group"),
                             br(),
                             br(),
                             downloadButton("downloaddel", "Download csv file"),
                             br(),
                             br(),
                             dataTableOutput('deletion')
                        )
                        
               )
    )
  )
  
  server <- function(input, output) {
    #table
    output$table_compare <- renderDataTable(uncommun_simple)
    # Downloadable csv
    
    output$downloadData <- downloadHandler(
      filename = function() {
        paste("snp_patient_control", '.csv', sep='')
      },
      content = function(file) {
        write.csv(uncommun_simple, file, row.names = FALSE)
      }
    )
    #patient
    output$TiTv_patient <- renderPlot({
      counts_titv <- table (uncommun$TiTv_patient)
      ggplot(as.data.frame(table(uncommun$TiTv_patient)), aes(x=Var1,y=Freq,fill=Var1))+
        geom_bar(stat = 'identity')+labs(x="",y="Mutations",fill="")+
        theme(legend.position = "none",
              text = element_text(size = 20))+
        ylab("SNP count")
    })
    #count snptype
    output$snp_patient <- renderPlot({
      counts_snp <- table (uncommun$snp_patient)
      ggplot(as.data.frame(table(uncommun$snp_patient)), aes(x=Var1,y=Freq,fill=Var1))+
        geom_bar(stat = 'identity')+labs(x="",y="Mutations",fill="")+
        theme(legend.position = "none",
              text = element_text(size = 20))+
        ylab("SNP count")
    })
    #
    output$insertion <- renderDataTable(addition)
    # Downloadable csv
    
    output$downloadins <- downloadHandler(
      filename = function() {
        paste("insertion", '.csv', sep='')
      },
      content = function(file) {
        write.csv(addition, file, row.names = FALSE)
      }
    )
    #
    output$deletion <- renderDataTable(deletion)
    # Downloadable csv
    
    output$downloaddel <- downloadHandler(
      filename = function() {
        paste("deletion", '.csv', sep='')
      },
      content = function(file) {
        write.csv(deletion, file, row.names = FALSE)
      }
    )
  }
  
  # Run the application
  shinyApp(ui = ui, server = server)
}

