.onAttach <- function() {
  packageStartupMessage(paste0("This package was created with the Linux operating system in mind.", 
                        "Functions callcnv(), callindels(), and callsnp() are nonfunctional on other operating systems."))
}
