# Pipeline validation Workflow 
To validate the pipelines and demonstrate the utility of the developed tool, we performed a case study on the SRP074425 dataset.
### Expression data
FASTQ files obtained from the SRA database were first preprocessed. Quality was controlled and then reads were aligned to the Homo sapiens genome, assembly hg19. 
The gene count was determined using the gene_counting() function and the output count csv file was visualized using the EXPviz() shiny app function. First the count data was visualized before the normalization and after normalization. 

After assessing the sample’s variability, we analyzed gene expression. The expression analysis is performed by the app and the genes expression profile is represneted in pca plot, volcano plot, ma plot, and heatmap.
The  PCA of the dataset showed that the samples were clustered according to the disease (BM vs control based on the second component PC2).

When comparing control and adenocarcinoma brain metastasis (BM) samples, we detected 1757 differentially expressed genes (DEG), of which 718 were upregulated and 1039 were downregulated genes (p-value < 0.05 and logFC >2).
The expression is represented in the MA,volcano plots, and heatmaps.

### Ontology
The ontology analysis is performed by the EXPviz app and the top 20 ontologies (BP, CC, MF) and  are represented in Barplot and dotplot.


## SNPs Calling
The BAM files generated from the fastqProcession() function are used to call the SNPs with call_snp() function and the output is a VCF file containing the SNPs data. The output was visualized using the SNPviz() shiny app function.

SNP types are first determined in comparison with the reference genome and then are represented in a pie plot. Transition (Ti) and Transversion (Tv) distribution are represented in a barplot.

The variants in the coding region are calculated, then the nonsense, nonsynonymous, and synonymous mutation types are determined and represented in a barplot.

## CNVs Calling
The call_cnv() function was used to call CNVs from the preprocessed BAM files. The output is a csv file containing the copy number states of genomic regions.
Recurrent CNV regions are identified and represented in barplots. Then the original CNV calls on overlapping genomic features (protein-coding genes) are illustareted is an oncoPrint plot.The overlap permutation test was performed and presented in a plot.

## INDELs Calling
The call_indel() was used to call INDELs from the preprocessed BAM files. The output is a VCF file containing the INDELs data.

