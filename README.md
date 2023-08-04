# Exvar: A gene expression and genetic variation data analysis and visualization R package

## Table of Contents
1. [Motivation](#Motivation)
2. [About Exvar](#About-Exvar)
3. [Supported Species](#Supported-Species)
4. [Installation](#Installation)
5. [Functions](#Functions)
6. [License](#License)
7. [Contact](#Contact)
8. [Citation](#Citation)
9. [Contributors](#Contributors)


## Motivation

RNA sequencing data manipulation workflows are complex and require various skills and tools. This creates the need for user-friendly and integrated genomic data analysis and visualization tools.

We developed a novel R package using multiple R packages to perform gene expression analysis and genetic variant calling from RNA sequencing data. The package could be used to analyze eight species’ data. Multiple public datasets were analyzed using the developed package to validate the pipeline for all the supported species.

## About Exvar

The Exvar R package performs gene expression and  genetic variats (SNPs, Indels, and CNVs) data analysis and  visualization.

It is developed to facilitate and minimize the skills required for the genetic expression and variants calling. It includes 9 functions, providing an easy biologist-friendly workflow.

## Supported Species

It could be used to analyze data from eight species including:

    * Homo Sapiens (19, 37, and 38)

    * Mus Musculus

    * Arabidopsis Thaliana

    * Drosophila Melanogaster

    * Danio rerio

    * Rattus norvegicus

    * Caenorhabditis elegans  

    * and Saccharomyces Cerevisiae.

## Installation

The package could be installed as follows:

    install.packages("devtools")

   lLibrary(devtools)

    devtools::install_github("omicscodeathon/Exvar/Package")

   lLibrary(Exvar)

## Functions

    - requirement() >> Install required packages

    - processfastq() >>Preprocess fastq files

    - counts() >> gene Count analysis

    - expression() >>Identify DEGs

    - callsnp() >> SNP calling

    - callcnv() >> CNV calling

    - callindel() >> Indel calling

    - vizexp() >> Analyze and visualize gene expression data

    - vizsnp() >> Analyze and visualize SNP data

    - vizcnv() >> Analyze and visualize CNV data


## License  

Artistic license 2.0


## Citation

Hiba Ben Aribi, Imraan Dixon , Najla Abassi, and  Olaitan I. Awe. Exvar: An R Package for Gene Expression And Genetic Variation Data Analysis And Visualization. (2023) https://github.com/omicscodeathon/Exvar

## Contributors

   - Hiba Ben Aribi, Faculty of Sciences of Tunis, University of Tunis El Manar, Tunis, Tunisia. [![Linkdin URL](https://www.linkedin.com/in/hiba-ben-aribi-b29634239/)

   - Imraan Dixon, Faculty of Health Sciences, University of Cape Town, Cape Town, South Africa.

   - Najla Abassi, Higher Institute of Biotechnology Sidi Thabet, Manouba University, Tunisia.

   - Olaitan I. Awe, African Society for Bioinformatics and Computational Biology, Cape Town, South Africa.
