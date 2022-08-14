# Raw data Analysis Pipelines

## Fastq files preprocessing

- Read fastq files using "ShortRead" package
- Quality Control
    - Filter low quality reads using the "Rfastp" package
    - Extract QC data using "qcsummary()" function
    - Plot QC aspect using the "CurvePlot(xxx.json)" function
- Alignment 
    - Retreive reference genome information from the " BSgenome" library
    - Creat a DNAStringSet object
    - Creat FASTA file using "writeXStringSet()"
    - Build index using the "Buildindex()" function
    - alignment using the "Rsubread" package - "subjunc" algorithm
    [ require an annotation SAF file; retreive exon location using the "exons()" function on the refernce genome then coverted to a dataframe using "data.frame()" function]
    - mapping using the "subjunc()" function 
- Sort and index reads using "RsamTools" packages "sortBam()" and "indexBam()" functions to sort and index the files repectively
- read in IGV
- gene Counting
    - create GRange object using the "exonsby()"  function
    - read bam files using the "BamFile()" function
    - count reads in bam files using the "summarizeOverlaps()" function
    - retrieve a matrix of counts from either RangedSummarizedExperiment object using the assay() function
    
## Gene expression
- using "Rsamtools" package
- Organize bam files using the "BamfileList()" function
- creat "RangedSummarizedExperiment" object 
    - "exonsBy()"
    - "summarizeOverlaps()"
- parallelization. using the "BiocParallel" package
- "assay()" function
- rowRanges() 

## Differentially expressed genes
-  "DESeq2" package
    - create dataframe using the "data.frame()" function
    - create DESeq2 object
        - assay()
        - rowRanges()
        - DESeqDataSetFromMatrix()
        - create RangedSummarizedExperiment object 
        - DESeqDataSet() function 
    - DESeq() function
    - retrieve normalized and unnormalized values from our DESeq2 object using the counts() function [normalized = TRUE]
    - review the Variance/mean relationship and shrinkage using plotDispEsts()
    - extract our contrast of interest using the results() function
    - sort by p value DESeqResults[order(myRes$pvalue), ]
    - summary()
    - review the relationship between fold-changes and expression levels with a MA-Plot by using the plotMA() function and DESeqResults object
    - downweight genes with high fold change but low significance using the lfcShrink() function
    - convert the DESeqResults objects into a standard data frame using the as.data.frame function.
    
### Multiple testing
we use the Benjamini-Hockberg = ((p-value of gene)*(Total genes tested))/(rank of p-value of gene)
- add adjusted p-values column using the p.adjust() function
- filter genes with NA padj values using  myResAsDF[!is.na(myResAsDF$padj), ]


###Annotation
- retrieve Gene Symbols for the Entrez IDs using the select() function from the "org.db" packages
- merge the Entrez ID to Symbol table into  differential expression results table using the "merge()" function

## Genomic Variants (CNVs)

    



    



    
    
