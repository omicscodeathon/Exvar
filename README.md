# Exvar: A Gene Expression And Genetic Variants Data Analysis And Visualization R package

## Introduction

The rapid development of high-throughput sequencing technology represents great opportunities for researchers to investigate genomic variations and transcriptomic changes. Multiple sequencing data are publicly available on multiple databases. But unfortunately, the analysis techniques are complex and require multiple skills.


## Solution = Exvar R package

The Exvar R package performs gene expression and  genetic variats (SNPs, Indels, and CNVs) data analysis and  visualization.

It is developed to facilitate and minimize the skills required for the genetic expression and variants calling. It includes 9 functions, providing an easy biologist-friendly workflow.

The package functions:

![image](https://user-images.githubusercontent.com/73958439/216056825-f40b8737-aa2e-41ee-ba70-540204ba1b6f.png)
![image](https://user-images.githubusercontent.com/73958439/216056956-3e7a2d7c-e126-4564-b638-70201b75275c.png)


Documentation: omicscodeathon/Exvar/Exvar.pdf

Examples of each function output: https://github.com/omicscodeathon/Exvar/tree/main/Supplementary_materials

## Installation

The package could be installed as follows:

install.packages("devtools")

devtools::install_github("omicscodeathon/Exvar/Package")

Library(Exvar)

## Workflow

 STEP1: Install the requirements using the requirement() function.

 STEP2: preprocess the fastq files using the processfastq() function.

 STEP3 : perform the gene counting using the counts() function and visualize the output csv file using the vizexp() function.

 STEP4: Call variants using the callSNP(), callCNV(), and callIndel() functions and visualize the data using the vizsnp() and vizcnv() functions.

 STEP5: download the results.

## Contributors

- Hiba Ben Aribi, Faculty of Sciences of Tunis, University of Tunis El Manar, Tunis, Tunisia.

- Imraan Dixon, Faculty of Health Sciences, University of Cape Town, Cape Town, South Africa.

- Najla Abassi, Higher Institute of Biotechnology Sidi Thabet, Manouba University, Tunisia.

- Olaitan I. Awe, African Society for Bioinformatics and Computational Biology, South Africa.
