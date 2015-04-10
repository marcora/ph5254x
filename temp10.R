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
