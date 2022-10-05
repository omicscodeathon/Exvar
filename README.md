# Gene Expression And  Genetic Variants Data Analysis And Visualization

![image](https://user-images.githubusercontent.com/73958439/193952491-dba21569-10ea-4e9e-b651-66f6cde66aa9.png)

## Introduction
The rapid development of high-throughput sequencing technology represents great opportunities for researchers to investigate genomic variations and transcriptomic changes. Multiple sequencing data are publicly available on multiple databases. But unfortunately, the analysis techniques are complexe and require multiple skills. 

## The project :

The project aims to create a workflow to perform gene expression and  genetic variats (SNPs, Indels, and CNVs) data analysis and  visualization.
The project includes a data analysis R package  and data visualization shiny app .

## 1. Data analysis R package : ExpVar

The NamePackage  R package is developed to facilitate the data analysis and minimize the skills required for the genetic expression and  variants calling. It only includes 5 functions, providing an easy biologist-friendly workflow. Standard analysis and  filtering settings are used as default

![r worklow](https://user-images.githubusercontent.com/73958439/193316510-27dbd891-a406-4f7f-a4a8-03c2c11ddbf2.png)


The package could be installed in R studio using this command :

install.packages("devtools")

devtools::install_github("omicscodeathon/EveQ/install_package")
The guidelines for implementing this tool are available in the dowument named "Package_user_guide.pdf"

The package consists of 6 Functions : preprocessing(); gene_counting() ; snps_calling(); cnvs_calling(); and indels_calling(). The package workflow is represented in Table1.


| Function  | Role| Input | Output | 
| ------------- | ------------- | ------------- | ------------- |
| preprocessing() |  preprocess raw data |  Fastq   | Bam  |
| gene_counting()  |   genes expression analysis  | Bam   | csv  |
| snps_calling()  |   SNPs calling  | Bam   | vcf  |
| indels_calling()  |   Indels calling  | Bam   | vcf  |
| cnvs_calling()  |   CNVs calling  | Bam   | vcf  |


## 2. Data visualization Shiny Apps

Three shiny apps were developed: EXPviz for genetic expression data visualization, SNPviz for SNPs data visualization, and CNVviz for CNVs  data visualization.

The apps integrate multiple analysis, however, with a biologist friendly web interface. They provides the possiblity to perform advance analysis and setting customization options without the need for any programming skills.

## 2.1. EXPviz

The EXPviz requires as input the genes count data csv file and a metadata file (defining the study design)  for the expression data visualization.  

The app includes 3 parts : 

(I) Count data visualization: the count data is visualized before and after the log transformation for quality control purposes. The genes count data are visualized in a density plot.  

(II) the expression analysis is performed, and the differentially expressed genes are defined. The user can customize the p-value and the LogFC value to be used. As default, the p-value is set to < 0.05 and the LogFC value to |logFC|>2 to identify the DEGs. 
The DEGs are then visualized in a volcano plot, an MA plot, a heatmap, and a PCA plot. The non-significant genes, the upregulated, and the downregulated gene are marked in different colors in the volcano plot. 

(III) the ontology analysis is performed, and the adjusted p-value<0.05 is set as default to identify significant ontologies. The p-value could be customized by the user.

The significant ontologies will be represented in tables and the top ontologies are represented in six barplots and dotplots: (i) BP associated with upregulated genes, (ii) MF associated with upregulated genes, (iii) CC associated with upregulated genes, (iv) BP associated with downregulated genes, (v) MF associated with downregulated genes, and (vi) CC associated with downregulated genes.

## 2.2. SNPviz

The app includes 3 parts : 

(I)  The SNPs types are determined in comparison with the reference genome (AtoC, AtoG, AtoT, CtoA, CtoG, CtoT, TtoC, TtoG, TtoA, GtoA, GtoC, GtoT), and the are represented in a pie chart.

(II) The SNPs are then divided into The Transition (Ti): purine-to-purine, pyrimidine-to-pyrimidine, and Transveersion (Tv): purine-to-pyrimidine, pyrimidine-to-purine. And the result is represented in a barplot.

(III) The Trinucleotide motif analysis results are represented in a barplot showing the trinucleotide distributions.

(III) amino Acids changes are determined and counted.

## 2.3. CNVviz
The app includes 3 parts : 

(I) Recurrent CNV regions
<p align="center">
  <![image](https://user-images.githubusercontent.com/73958439/193952944-ed598959-462a-4845-9571-7e1ea40b5b12.png) width="300"/>
</p>




(II) CNVs - genes overlapping

(III)overlap mutated test

## Installation

The ExpVar package could be installed in R studio using this command :

install.packages("devtools")

devtools::install_github("omicscodeathon/EveQ/ExpVar")

The  Web shiny apps  is publicly available at www.shinyapps.io/......

The guidelines for implementing these tool are available in the dowument named "Shinny_app_user_guide.pdf"

## Case study
To validate the pipelines and demonstrate the utility of the developed tool, we performed a case study investigating the Blood-Brain-barrier dysfunction mechanisms in the brain tumor patients.

The analyzed data and results are available in "Case_study_results.pdf"


## Acknowledgements
This project was performed during the October 2022 Omics codeathon organized in collaboration with the African Society for Bioinformatics and Computational Biology (ASBCB). The authors thank the National Institutes of Health (NIH) Office of Data Science Strategy (ODSS) and the National Center for Biotechnology Information (NCBI) for their immense support before and during the codeathon.

## Contributors

- Hiba Ben Aribi, Faculty of Sciences of Tunis, University of Tunis El Manar (UTM), Tunis, Tunisia

- Imraan Dixon, Faculty of Health Sciences, University of Cape Town (UCT), Cape Town, South Africa 

- Najla Abassi, Higher Institute of Biotechnology Sidi Thabet, Manouba University (UMA), Tunisia
