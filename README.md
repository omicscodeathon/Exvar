# Gene Expression And  Genetic Variants Data Analysis And Visualization

## Introduction
The rapid development of high-throughput sequencing technology represents great opportunities for researchers to investigate genomic variations and transcriptomic changes. Multiple sequencing data are publicly available on multiple databases. But unfortunately, the analysis techniques are complexe and require multiple skills. 

## The project : Shiny App for gene expression and genetic variants  data visualization

The project aims to create a workflow to perform gene expression and  genetic variats (SNPs, Indels, and CNVs) data analysis and  visualization.
The project includes a data analysis R package  and data visualization shiny app .

## The  data analysis R package 

The NamePackage  R package is developed to facilitate the data analysis and minimize the skills required for the genetic expression and  variants calling. It only includes 5 functions, providing an easy biologist-friendly workflow. Standard analysis and  filtering settings are used as default

![r worklow](https://user-images.githubusercontent.com/73958439/193316510-27dbd891-a406-4f7f-a4a8-03c2c11ddbf2.png)


The package could be installed in R studio using this command :

install.packages("devtools")

devtools::install_github("omicscodeathon/EveQ/install_package")

The package consists of 6 Functions : preprocessing(); gene_counting() ; snps_calling(); cnvs_calling(); and indels_calling(). The package workflow is represented in Table1.


| Function  | Role| Input | Output | 
| ------------- | ------------- | ------------- | ------------- |
| preprocess() |  preprocess raw data |  Fastq   | Bam  |
| expression()  |   genes expression analysis  | Bam   | csv  |
| SNPs()  |   SNPs calling  | Bam   | vcf  |
| Indels()  |   Indels calling  | Bam   | vcf  |
| CNVs()  |   CNVs calling  | Bam   | vcf  |


## VEveQ Shiny App for data visualization

The VEveQ RShiny package can be downloaded from  https://github.com/omicscodeathon/EVeQ/VEveQ_package. 

The VEVeQ Web shiny app  is publicly available at www.shinyapps.io/VECeQ.

The guidelines for implementing this tool and related updates, are available at:  https://github.com/omicscodeathon/EVeQ/VEVeQ_Shinyapp/README.md.

## Case study
To demonstrate the utility of the developed tool, we performed a case study investigating the different Blood-Brain-barrier dysfunction mechanisms, in different neurological diseases Huntington’s disease (HD),  neuroinfectious disease and the brain tumor pathology.

The analyzed data and results are available in : https://github.com/omicscodeathon/ECeQ/case_study.

## Abbreviations
- EVeQ : expression, genetic Variants and eQTL analysis 

- VEVeQ : Visualize expression, CNV and eQTL data

- DEGs : Differentially Expressed Genes

- CNV : Copy Number Variation

- eQTL : Expressed Quantitative Trait Loci

## Acknowledgements
This project was performed during the October 2022 Omics codeathon organized in collaboration with the African Society for Bioinformatics and Computational Biology (ASBCB). The authors thank the National Institutes of Health (NIH) Office of Data Science Strategy (ODSS) and the National Center for Biotechnology Information (NCBI) for their immense support before and during the codeathon.

## Contributors

- Hiba Ben Aribi, Faculty of Sciences of Tunis, University of Tunis El Manar (UTM), Tunis, Tunisia

- Imraan Dixon, Faculty of Health Sciences, University of Cape Town (UCT), Cape Town, South Africa 

- Najla Abassi, Higher Institute of Biotechnology Sidi Thabet, Manouba University (UMA), Tunisia
