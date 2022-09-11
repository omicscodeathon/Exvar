#CNV calling
print("Calling copy number variants...")
variant_df <- 
  vcf2txt(paste0(file_path_sans_ext(basename(inputFastq[c(x)])), ".vcf"))
cnv_data <- cnv.data(variant_df)
cnv_call <- cnv.call(cnv_data)
cnvcf <- asVCF(cnv_call, 
               paste0(file_path_sans_ext(basename(inputFastq[c(x)]))))
writevcf(cnvcf, paste0(file_path_sans_ext(basename(inputFastq[c(x)])), 
                       ".vcf"),
         index = TRUE)
