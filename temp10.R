<<<<<<< HEAD
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
=======
# library(devtools)
# install_github("genomicsclass/ERBS")
library(ERBS)

data(HepG2)
data(GM12878)

fo = findOverlaps(HepG2, GM12878)
erbs = HepG2[queryHits(fo)]
erbs = granges(erbs)

unique(genome(erbs))

library(BSgenome.Hsapiens.UCSC.hg19)

erbsseq = getSeq(BSgenome.Hsapiens.UCSC.hg19, erbs)

gccontent = sapply(erbsseq, function(seq){
  
})

gccontent = letterFrequency(erbsseq, "GC") / width(erbsseq)
median(gccontent)

rerbsseq = getSeq(BSgenome.Hsapiens.UCSC.hg19, shift(erbs, 10000))
gccontent = letterFrequency(rerbsseq, "GC") / width(rerbsseq)
median(gccontent)
>>>>>>> f24ad0feac8ed22e87b68bbd07eec4c2e257ce55
