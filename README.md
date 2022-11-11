# Gene Expression And  Genetic Variants Data Analysis And Visualization

![image](https://user-images.githubusercontent.com/73958439/193952491-dba21569-10ea-4e9e-b651-66f6cde66aa9.png)

## Introduction
The rapid development of high-throughput sequencing technology represents great opportunities for researchers to investigate genomic variations and transcriptomic changes. Multiple sequencing data are publicly available on multiple databases. But unfortunately, the analysis techniques are complex and require multiple skills. 

## The project :

The project aims to create a workflow to perform gene expression and  genetic variats (SNPs, Indels, and CNVs) data analysis and  visualization.


## ExpVar R package :

The ExpVar  R package performs gene expression and  genetic variats (SNPs, Indels, and CNVs) data analysis and  visualization.
IT is developed to facilitate and minimize the skills required for the genetic expression and  variants calling. It includes 8 functions, providing an easy biologist-friendly workflow. 

The package functions are represented in Table1.

![image](https://user-images.githubusercontent.com/73958439/194703211-02d5b899-f51a-4c6d-906d-d4ada8d6b570.png)

The apps interfaces are documented in:

  * EXPviz: 
 
 https://github.com/omicscodeathon/ExpVar/blob/main/Supplementary_materials/EXPviz.pdf
 
  * SNPviz: 
 
 https://github.com/omicscodeathon/ExpVar/blob/main/Supplementary_materials/SNPviz.pdf
 
  * CNVviz: 
 
 https://github.com/omicscodeathon/ExpVar/blob/main/Supplementary_materials/CNVviz.pdf


## Install 

The package could be installed in R studio using this command :

install.packages("devtools")

devtools::install_github("omicscodeathon/ExpVar/ExpVar/ExpVar_package")

## Workflow

Note: gmapR::GmapGenome() replaces createReference() functionality
![image](https://user-images.githubusercontent.com/73958439/194729911-b0c74261-0a3f-4fcf-8608-454fec8592c6.png)



## Acknowledgements
This project was performed during the October 2022 Omics codeathon organized in collaboration with the African Society for Bioinformatics and Computational Biology (ASBCB). The authors thank the National Institutes of Health (NIH) Office of Data Science Strategy (ODSS) and the National Center for Biotechnology Information (NCBI) for their immense support before and during the codeathon.

## Contributors

- Hiba Ben Aribi, Faculty of Sciences of Tunis, University of Tunis El Manar (UTM), Tunis, Tunisia

- Imraan Dixon, Faculty of Health Sciences, University of Cape Town (UCT), Cape Town, South Africa 

- Najla Abassi, Higher Institute of Biotechnology Sidi Thabet, Manouba University (UMA), Tunisia
