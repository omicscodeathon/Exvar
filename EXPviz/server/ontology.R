

logFC_input = 0.05
p_value  =0.05
countdata <- as_data_frame(read_csv("demo/Expression_BT.csv")) 
expression = countdata
expression$diffexpressed <- "Not Significant"
expression$diffexpressed[expression$log2FoldChange  > logFC_input  & expression$pvalue < p_value ] <- "Upregulated genes"
expression$diffexpressed[expression$log2FoldChange  < -logFC_input  & expression$pvalue < p_value ] <- "Downregulated genes"


#filter upregulated genes
upgene <- as.data.frame(expression[expression$diffexpressed == "Upregulated genes",])   # output 
upgenes_list <- upgene$ENSEMBL

downgene <- as.data.frame(expression[expression$diffexpressed == "Downregulated genes",])   # output 
downgenes_list <- downgene$ENSEMBL

#BP_up   
#output$BP_up_table<- renderPlot({
BP_up_go <- enrichGO(gene = upgenes_list, OrgDb="org.Hs.eg.db", keyType="ENSEMBL", ont="BP")
BP_up_table <- as.data.frame(BP_up_go)
#  BP_up_table
#})

#input number ontologies in plot
#showcat <- input$show
#gnelist <- if  list("all DEGs" = 1, "Upregulated genes" = 2, "Downregulated genes" = 3) #output$value <- renderPrint({ input$select })
output$BP_barplot<- renderPlot({
  barplot_BP_up <-plot(barplot(BP_up_go, showCategory = 20))
  barplot_BP_up

})
output$BP_dotplot<- renderPlot({
  dotplot(BP_up_go, showCategory=30)
})

####
#output$MF_up_table<- renderPlot({
  MF_up_go <- enrichGO(gene = upgenes_list, OrgDb="org.Hs.eg.db", keyType="ENSEMBL", ont="MF")
#})
output$MF_barplot<- renderPlot({
  barplot_MF_up <-plot(barplot(MF_up_go, showCategory = 20))
  barplot_MF_up
})

output$MF_dotplot<- renderPlot({
  dotplot(BP_up_go, showCategory=30)
})

###
#output$MF_up_table<- renderPlot({
CC_up_go <- enrichGO(gene = upgenes_list, OrgDb="org.Hs.eg.db", keyType="ENSEMBL", ont="CC")
#})
output$MF_barplot<- renderPlot({
  barplot_CC_up <-plot(barplot(CC_up_go, showCategory = 20))
  barplot_CC_up
})

output$CC_dotplot<- renderPlot({
  dotplot(CC_up_go, showCategory=30)
})

##
