# Workflow

## Expression data
FASTQ files obtained from the SRA database were first preprocessed. Quality was controlled and then reads were aligned to the Homo sapiens genome, assembly hg19. 
The gene count was determined using the gene_counting() function and the output count csv file was visualized using the EXPviz() shiny app function. First the count data was visualized before the normalization and after normalization. 

After assessing the sample’s variability, we analyzed gene expression. The expression analysis is performed by the app and the genes expression profile is represneted in pca plot, volcano plot, ma plot, and heatmap.
The  PCA of the dataset showed that the samples were clustered according to the disease (BM vs control based on the second component PC2).

When comparing control and adenocarcinoma brain metastasis (BM) samples, we detected 1757 differentially expressed genes (DEG), of which 718 were upregulated and 1039 were downregulated genes (p-value < 0.05 and logFC >2).
The expression is represented in the MA,volcano plots, and heatmaps.

### Ontology
The ontology analysis is performed by the EXPviz app and the top 20 ontologies (BP, CC, MF) and  are represented in Barplot and dotplot.

Figures are in EXPviz.pdf

