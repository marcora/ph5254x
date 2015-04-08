library(Biobase)
sessionInfo() 

source("http://www.bioconductor.org/biocLite.R")
biocLite()

library(devtools)
install_github("genomicsclass/ph525x")

# BSgenome.Hsapiens.UCSC.hg19
# SNPlocs.Hsapiens.dbSNP.20120608
# Homo.sapiens

# library(BiocInstaller)
# biocLite([missing package name here, in quotes])
library(BiocInstaller)
biocLite("BSgenome.Hsapiens.UCSC.hg19")
