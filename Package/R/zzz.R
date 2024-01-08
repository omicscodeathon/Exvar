.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0("This package was created with the Linux operating system in mind.", 
                        "Functions callcnv(), callindels(), callsnp(), count(), and geneExpression() are nonfunctional on other operating systems.\n",
                              "Additionally, a samtools installation is recommended for snp and indel calling."))
}
