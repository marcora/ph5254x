# source("http://bioconductor.org/biocLite.R")
# biocLite("affy")

library(affy)

datadir <- "rawdata-master"
basedir <- paste0(datadir, "/celfiles")

tab = read.delim(file.path(basedir, "sampleinfo.txt"), check.names=FALSE, as.is=TRUE)
rownames(tab) = tab$filenames
tab

fns = list.celfiles(basedir)
fns

all(fns %in% tab[,1]) # check for consistency

ab = ReadAffy(filenames=file.path(basedir, rownames(tab)), phenoData=tab)

dim(pm(ab))
