# source("http://bioconductor.org/biocLite.R")
# biocLite("Gviz")

# library(devtools)
# install_github("genomicsclass/ph525x")
library(ph525x)

modPlot("ESR1", useGeneSym=FALSE, collapse=FALSE)

keytypes(Homo.sapiens)
columns(Homo.sapiens)
select(Homo.sapiens, keytype="SYMBOL", key="ESR1", columns="GENEID")

library(Homo.sapiens)
length(transcriptsBy(Homo.sapiens, by="gene")$"2099")

