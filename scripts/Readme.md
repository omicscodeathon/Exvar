# Raw data Analysis Pipelines
# Note : This is not the final pipeline !!
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
- parallelization using the "BiocParallel" package
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


### Annotation
- retrieve Gene Symbols for the Entrez IDs using the select() function from the "org.db" packages
- merge the Entrez ID to Symbol table into  differential expression results table using the "merge()" function

## Genomic Variants (CNVs)
- variant calling using the "Bcftools" package: bam to vcf file
- read vcf files using the readVcf() function
- annotation using the "VariantAnnotation" package 
- The variant information is recorded in a VRange object; rowRanges(vcf)
- retrive dbSNP data from refernce genome (SNPlocs.Hsapiens.dbSNP144.GRCh37)
        - Extract chromosome name from the VRange object: seqnames()@values
        - Retrive SNPs by chromosome: snpsBySeqname()
- update the seqlevelsStyle to "UCSC"

    #change seqlevelsStyle
    
    seqlevelsStyle(my_snps) <- "UCSC"
    
    #change genome
    
    genome(my_snps) <- "hg19"
    
- make a rsID table
- Extract information of variants using the gsub() function
- marge the variant and rsID data.frames using the "merge()" function
- predict amino acid changes using the predictcoding() function
- annotate and determnie mutation types
- integrate SNP and amino acid changes

## eQTL
- quality control
         - remove markers with excessive missing genotypes
         - exclude subjects with excessive missing genotypes
         - exclude subjects with gender-mismatch
         - remove markers with heterozygous haploid genotypes
         - remove markers violating Hardy-Weinberg Equilibrium (HWE)
         - remove markers with informative missingness
         - remove markers with low minor allele frequency (MAF)
         - exclude subjects with abnormal heterozygosity rate
         - exclude related subjects
         - identify and exclude individuals with divergent ancestry
- eQTL mapping using the " MatrixEQTL" package
- cis-eQTL association
- trans-eQTL association
    



    
    
