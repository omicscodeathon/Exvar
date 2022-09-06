# AECeQ

## About
###Fastq pre-processing
The package makes use of the “Rfastp” package to perform quality control analysis on one or more fastq files (Wang & Carroll, 2022) . It filters out bad reads and 
trims adaptors, producing a “quality_checked” fastq file and a quality control summary. The “gmapR” package is used to build an index from a reference genome of which 
the quality controlled fastq files are aligned to (Barr et al., 2022). Bundled within gmapR are functions relating to the Genomic Short-read Nucleotide Alignment 
Programme (GSNAP) tool. This allows alignment of the fastq files to the indexed reference genome. The output is an indexed bam file.



## Install 

The package could be installed in R studio using this command :

devtools::install_github("omicscodeathon/ECeQ/AECeQ/AECeQ_package")

## Requirements
.
.
.
.

## Functions

| Function  | Role| Input | Output | 
| ------------- | ------------- | ------------- | ------------- |
| preprocess() |  preprocess raw data |  Fastq   | Bam  |
| expression()  |   genes expression analysis  | Bam   | csv  |
| variant()  |   variants calling  | Bam   | vcf  |
| eQTL()  |   eQTL analysis  | Bam   |     |
| ECeQ()  |   all analysis  | Fastq   | VECeQ interface  |


## Citation
.
.
.
