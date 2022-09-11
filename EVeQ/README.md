# EVeQ

EVeQ ia a  Gene Expression, Genomic Variants and eQTL Data  Analysis R package

## Install 

The package could be installed in R studio using this command :

devtools::install_github("omicscodeathon/ECeQ/AECeQ/AECeQ_package")

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


