# EVeQ

EVeQ ia a  Gene Expression, Genomic Variants and eQTL Data  Analysis R package

## Install 

The package could be installed in R studio using this command :

devtools::install_github("omicscodeathon/ECeQ/AECeQ/AECeQ_package")

## Workflow


## Functions

| Function  | Role| Input | Output | 
| ------------- | ------------- | ------------- | ------------- |
| preprocess() |  preprocess raw data |  Fastq   | Bam  |
| expression()  |   genes expression analysis  | Bam   | csv  |
| variant()  |   variants calling  | Bam   | vcf  |
| eQTL()  |   eQTL analysis  | Bam   |     |
| ECeQ()  |   all analysis  | Fastq   | VECeQ interface  |

