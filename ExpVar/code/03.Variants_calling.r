## Variant calling via VariantTools
print("Calling variants...")
parallelparam <- MulticoreParam(workers = 4)
genome(repeats) <- genome(refgen)
tallyParam <- TallyVariantsParam(refgen, mask = repeats, BPPARAM = parallelparam)
variantcall <- callVariants(bamfl, tallyParam)
print("Variant calling completed...")
sampleNames(variantcall) <- 
  paste0(file_path_sans_ext(basename(inputFastq[c(x)])))
mcols(variantcall) <- NULL
vcf <- asVCF(variantcall)
print("Writing vcf...")
writeVcf(vcf, paste0(file_path_sans_ext(basename(inputFastq[c(x)])), ".vcf"), 
         index = TRUE)
