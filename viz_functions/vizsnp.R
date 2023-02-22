#'  vizsnp
#'
#' @param dir Path contain two folders named "Control" and "patient" containing the VCF files.
#'
vizsnp <- function( dir) {
#load requirement

library(shiny)
library(ggplot2)
library(data.table)
library(fs)
library(shinythemes)
library(dplyr)
library(purrr)
library(shinydashboard)
library(DT)
###function input : dir
setwd(dir)
#file_paths_control
file_paths_control <- fs::dir_ls("control")
#file_paths_patient
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
comined_control <-file_contenets_control %>% map(as_tibble) %>%
  reduce(full_join)
# delete NAs
# result : keep commun between all samples
comined_control <- na.omit(comined_control)
## identify snps
# Transition (Ti)
ti <- c("A>G","G>A","C>T","T>C")
# Transveersion (Tv)
tv <- c("A>T","A>C","G>T","G>C","C>A","C>G","T>A","T>G")
#
comined_control$nuSub <- paste0(comined_control$REF,">",comined_control$ALT)
comined_control$TiTv[comined_control$nuSub %in% ti] <- "Ti"
comined_control$TiTv[comined_control$nuSub %in% tv] <- "Tv"
#
comined_control <-dplyr::rename(comined_control,"snp_controls"="nuSub")
comined_control <-dplyr::rename(comined_control,"TiTv_controls"="TiTv")
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
comined_patient$nuSub <- paste0(comined_patient$REF,">",comined_patient$ALT)
comined_patient$TiTv[comined_patient$nuSub %in% ti] <- "Ti"
comined_patient$TiTv[comined_patient$nuSub %in% tv] <- "Tv"
#
comined_patient<-dplyr::rename(comined_patient,"snp_patient"="nuSub")
comined_patient<-dplyr::rename(comined_patient,"TiTv_patient"="TiTv")

# Part III: compare groups
compare_group <- merge(comined_patient, comined_control, by = c("#CHROM", "POS"),all = TRUE)
compare_group[is.na(compare_group)] <- "no"
compare_group$compare[compare_group$snp_controls !="no" &  compare_group$snp_patient =="no" ] <- "deletion"
compare_group$compare[compare_group$snp_controls =="no" &  compare_group$snp_patient !="no" ] <- "addition"
compare_group$compare[compare_group$snp_controls ==  compare_group$snp_patient ] <- "same"
compare_group$compare[compare_group$snp_controls !=  compare_group$snp_patient &  compare_group$snp_controls !="no" &compare_group$snp_patient !="no"] <- "different"
uncommun <- compare_group[compare_group$compare == "different",]
commun <- compare_group[compare_group$compare == "same",]
#
addition <- compare_group[compare_group$compare == "addition",]
addition <- addition[c("#CHROM","POS","REF.x","ALT.x","TiTv_patient")]
addition<-dplyr::rename(addition,"Chromosome"="#CHROM")
addition<-dplyr::rename(addition,"TiTv"="TiTv_patient")
addition<-dplyr::rename(addition,"SNP position"="POS")
addition<-dplyr::rename(addition,"REF"="REF.x")
addition<-dplyr::rename(addition,"ALT"="ALT.x")
deletion <- compare_group[compare_group$compare == "deletion",]
deletion <-deletion[c("#CHROM","POS","REF.y","ALT.y","TiTv_controls")]
deletion<-dplyr::rename(deletion,"Chromosome"="#CHROM")
deletion<-dplyr::rename(deletion,"TiTv"="TiTv_controls")
deletion<-dplyr::rename(deletion,"SNP position"="POS")
deletion<-dplyr::rename(deletion,"REF"="REF.y")
deletion<-dplyr::rename(deletion,"ALT"="ALT.y")
# dataframe to shown
uncommun_simple <- data.frame(uncommun$`#CHROM`, uncommun$POS, uncommun$ID.x,
                              uncommun$REF.x, uncommun$ALT.x, uncommun$ALT.y)
#rename columns
names(uncommun_simple) <- c('CHROM','POS','ID','REF','ALT_patient', 'ALT_control')
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
#
##write.vcf(comined_patient, file="test.vcf")

