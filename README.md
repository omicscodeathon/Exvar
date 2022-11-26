# Gene Expression And  Genetic Variants Data Analysis And Visualization

## Introduction
The rapid development of high-throughput sequencing technology represents great opportunities for researchers to investigate genomic variations and transcriptomic changes. Multiple sequencing data are publicly available on multiple databases. But unfortunately, the analysis techniques are complex and require multiple skills. 


## Solution = ExpVar R package :

The ExpVar  R package performs gene expression and  genetic variats (SNPs, Indels, and CNVs) data analysis and  visualization.

It is developed to facilitate and minimize the skills required for the genetic expression and variants calling. It includes 9 functions, providing an easy biologist-friendly workflow. 

The package functions :

![image](https://user-images.githubusercontent.com/73958439/194703211-02d5b899-f51a-4c6d-906d-d4ada8d6b570.png)

Documentation : https://github.com/omicscodeathon/ExpVar/blob/main/ExpVar/Function%20Documentation.pdf

Examples of each function output : https://github.com/omicscodeathon/ExpVar/tree/main/Supplementary_materials

## Installation

The package could be installed as follows:

install.packages("devtools")

devtools::install_github("omicscodeathon/ExpVar/ExpVar/ExpVar_package")

Library(ExpVar)

## Workflow

 STEP1: preprocess the fastq files using the Preprocession() function.
 
 STEP2 : perform the gene counting using the gene_counting() function and visualize the output csv file using the EXPviz() function.
 
 STEP3: Call variants using the callSNP(), callCNV(), and callIndel() functions and visualize the data using the SNPviz() and CNVviz() functions.
 
 STEP4: download the results.


## Acknowledgements

This project was performed during the October 2022 Omics codeathon organized in collaboration with the African Society for Bioinformatics and Computational Biology (ASBCB). The authors thank the National Institutes of Health (NIH) Office of Data Science Strategy (ODSS) and the National Center for Biotechnology Information (NCBI) for their immense support before and during the codeathon.

## Contributors

- Hiba Ben Aribi, Faculty of Sciences of Tunis, University of Tunis El Manar, Tunis, Tunisia.

- Imraan Dixon, Faculty of Health Sciences, University of Cape Town, Cape Town, South Africa.

- Najla Abassi, Higher Institute of Biotechnology Sidi Thabet, Manouba University, Tunisia.
