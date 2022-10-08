# read demo  files  // if demo = true
metadata <- as.data.frame(read_table("demo/metadata.txt"))  
countdata <- as_data_frame(read_csv("demo/count.csv"))  

# read input files
#metadata <-  input$  # change to input file
#countdata <- input$  # change to input file

coldata <- metadata
countData = as.data.frame(countdata)  # change to input
rownames(countData) = countData$gene   #first column name
countData = countData[,-c(1)]
counts_data <- countData
# creat deseq dataset object  ## add button Run DE
factor = ~ factor1  #input
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = coldata,
                              design = factor)  # change factor to input  // + possibilty to add multiple

dds


# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10  # reads number could be customized
dds <- dds[keep,]

dds

# Step 3: Run DESeq ----------------------
dds <- DESeq(dds)
expression_data  <- results(dds)
#summary(expression_data )

# Explore Results ----------------
#p_value <- 0.01 # alpha = input
#expression_data0.01 <- results(dds, alpha = p_value)  
#summary(expression_data0.01 )   # show as output

# contrasts 
#resultsNames(dds)

######## ------------- PLOTS ---------------############
# MA plot  ==> change to 0.05   == input
output$ma_plot <- renderPlot({
  DESeq2::plotMA(expression_data0.01,  ylim = c(-5, 5))
})

# identify DEGss
#use p value or adjust # if loop #input
#demo
p_value = 0.01 
logFC_input = 2 
#input
#p_value <- input$pval
#logFC_input <- input$log
# add a column 
expression = as.data.frame(expression_data)
expression$diffexpressed <- "Not Significant"
expression$diffexpressed[expression$log2FoldChange  > logFC_input  & expression$pvalue < p_value ] <- "Upregulated genes"
expression$diffexpressed[expression$log2FoldChange  < -logFC_input  & expression$pvalue < p_value ] <- "Downregulated genes"

# volcano plot
output$volcano_plot <- renderPlot({
  ggplot(data=expression, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed)) + geom_point() + theme_minimal()
  
})

#pca
output$pca_plot <- renderPlot({
  rld = rlog(dds)
  factorpca <- "factor1"  #input
  plotPCA(rld, intgroup=factorpca )
  
})
# heatmap all genes
#p_value = 0.01 #input
#logFC_input = 2 #☻input
output$heatmap <- renderPlot({
  mat <-counts(dds, normalized=T)[rownames(expression),]
  mat.z <- t(apply(mat,1,scale))
  cold = coldata
  rownames(cold) <- cold[,1]
  colnames(mat.z)<-rownames(cold)
  mat.z
  
  Heatmap(mat.z, cluster_rows = T, cluster_columns = T, column_labels = colnames(mat.z), name="Z-score",
  ) # add row_labels = expression[rownames(mat.z),]$symbol  to get genes in map
  
})
# heatmap only DEGs
output$heatmap_degs <- renderPlot({
  #filter only DEGs 
  degs <- expression[(expression$pvalue>p_value)&(abs(expression$log2FoldChange)>logFC_input),]
  # heatmap all genes
  mat <-counts(dds, normalized=T)[rownames(degs),]
  mat.z <- t(apply(mat,1,scale))
  cold = coldata
  rownames(cold) <- cold[,1]
  colnames(mat.z)<-rownames(cold)
  mat.z
  
  Heatmap(mat.z, cluster_rows = T, cluster_columns = T, column_labels = colnames(mat.z), name="Z-score",
  ) # add row_labels = expression[rownames(mat.z),]$symbol  to get genes in map
  
})




