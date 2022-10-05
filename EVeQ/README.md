# ExpVar

ExpVar ia a  Gene Expression and Genomic Variant Analysis and Visualisation R package

## Install 

The package could be installed in R studio using this command :

install.packages("devtools")

devtools::install_github("omicscodeathon/ExpVar/ExpVar/ExpVar_package")

## Workflow
![image](https://user-images.githubusercontent.com/73958439/189546512-fa42662b-9005-47f4-94ac-91a11ef3b328.png)

## Functions

| Function  | Role| Input | Output | 
| ------------- | ------------- | ------------- | ------------- |
| preprocess() |  preprocess raw data |  Fastq   | Bam  |
| expression()  |   genes expression analysis  | Bam   | csv  |
| SNPs()  |   SNPs calling  | Bam   | vcf  |
| Indels()  |   Indels calling  | Bam   | vcf  |
| CNVs()  |   CNVs calling  | Bam   | vcf  |
| eQTL()  |   eQTL analysis  | Bam   |  vcf + csv |


