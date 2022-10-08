# load demo file
countData <- as_data_frame(read_csv("demo/count.csv")) 

#read input
#countData <- input$csv_file

# genes to row names
countdata = countData
rownames(countdata) = countData$gene   #first column name
countData = countData[,-c(1)]
# Raw data 

output$read_count1 <- renderPlot({
  par(mar=c(8,4,4,1)+0.1)
  barplot( colSums(countData)/1e6, las=3, main="Total read counts (millions)", ylab="Total read counts in millions")
  
})
## error
output$read_count2 <- renderPlot({
  hist(countData[,1], br=200, xlab="Number of Reads Counts per Feature", main="Histogram of Read Counts")
  
})
# log transformation + log transformed read count plot
logCountData = log2(1+countData)
par(mfrow = c(1, 2), mar=c(8,4,4,1))  # two columns

# UI ----
output$log_transformation1 <- renderPlot({
  hist(logCountData[,1], main="Histogram of Log Read Counts", xlab="Log transformed counts")
})
output$log_transformation2 <- renderPlot({
  boxplot(logCountData,las=3, main="Boxplot of Log Read Counts")
})
# density plot
output$density_plot <- renderPlot({
  logCountData = log2(1+countData)
  par(mfrow = c(1, 2), mar=c(8,4,4,1))  # two columns   change margins
  x <- logCountData
  myColors = rainbow(dim(x)[2])
  plot(density(x[,1]),col = myColors[1], lwd=2,
       xlab="Expresson values", ylab="Density", main= "Distribution of transformed data",
       ylim=c(0, max(density(x[,1])$y)+.02 ) )
  for( i in 2:dim(x)[2] )
    lines(density(x[,i]),col=myColors[i], lwd=2)
  legend("topright", cex=1.1,colnames(x), lty=rep(1,dim(x)[2]), col=myColors )	
})

#  scatter plot two samples ### no output ui yet ----------------------------------
#demo
sample1 = "SRR11426825"  
sample2 = "SRR11426823"  

#input
#sample1 <- input$sam1 
#sample2 <- input$sam2

output$scatter_plot <- renderPlot({
  plot(logCountData[,1],logCountData[,2], xlab=sample1, ylab=sample2) 
})
