# Gene Expression And  Genetic Variants Data Analysis And Visualization

## Introduction
The rapid development of high-throughput sequencing technology represents great opportunities for researchers to investigate genomic variations and transcriptomic changes. Multiple sequencing data are publicly available on multiple databases. But unfortunately, the analysis techniques are complexe and require multiple skills. 

## The project : Shiny App for gene expression and genetic variants  data visualization

The project aims to create a workflow to perform gene expression and  genetic variats (SNPs, Indels, and CNVs) data analysis and  visualization.
The project includes a data analysis R package  and data visualization shiny app .

## 1. Data analysis R package 

The NamePackage  R package is developed to facilitate the data analysis and minimize the skills required for the genetic expression and  variants calling. It only includes 5 functions, providing an easy biologist-friendly workflow. Standard analysis and  filtering settings are used as default

![r worklow](https://user-images.githubusercontent.com/73958439/193316510-27dbd891-a406-4f7f-a4a8-03c2c11ddbf2.png)


The package could be installed in R studio using this command :

install.packages("devtools")

devtools::install_github("omicscodeathon/EveQ/install_package")
The guidelines for implementing this tool are available in the dowument named "Package_user_guide.pdf"

The package consists of 6 Functions : preprocessing(); gene_counting() ; snps_calling(); cnvs_calling(); and indels_calling(). The package workflow is represented in Table1.


| Function  | Role| Input | Output | 
| ------------- | ------------- | ------------- | ------------- |
| preprocess() |  preprocess raw data |  Fastq   | Bam  |
| expression()  |   genes expression analysis  | Bam   | csv  |
| SNPs()  |   SNPs calling  | Bam   | vcf  |
| Indels()  |   Indels calling  | Bam   | vcf  |
| CNVs()  |   CNVs calling  | Bam   | vcf  |


## 2. Data visualization Shiny App 

We developed an interactive shiny app. The app integrate multiple analysis, however, it shows a biologist friendly web interface. It provides the possiblity to perform advance analysis and setting customization options without the need for any programming skills.

The Shiny package could be installed in R studio using this command :

install.packages("devtools")

devtools::install_github("omicscodeathon/EveQ/install_app")

The  Web shiny app  is publicly available at www.shinyapps.io/app

The guidelines for implementing this tool are available in the dowument named "Shinny_app_user_guide.pdf"

## Case study
To validate the pipelines and demonstrate the utility of the developed tool, we performed a case study investigating the Blood-Brain-barrier dysfunction mechanisms in the brain tumor patients.

The analyzed data and results are available in "Case_study_results.pdf"


## Acknowledgements
This project was performed during the October 2022 Omics codeathon organized in collaboration with the African Society for Bioinformatics and Computational Biology (ASBCB). The authors thank the National Institutes of Health (NIH) Office of Data Science Strategy (ODSS) and the National Center for Biotechnology Information (NCBI) for their immense support before and during the codeathon.

## Contributors

- Hiba Ben Aribi, Faculty of Sciences of Tunis, University of Tunis El Manar (UTM), Tunis, Tunisia

- Imraan Dixon, Faculty of Health Sciences, University of Cape Town (UCT), Cape Town, South Africa 

- Najla Abassi, Higher Institute of Biotechnology Sidi Thabet, Manouba University (UMA), Tunisia
