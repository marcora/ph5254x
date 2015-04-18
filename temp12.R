library(ERBS)

data(HepG2)

HepG2.hg19 = HepG2

library(rtracklayer)

ch = import.chain("hg19ToHg38.over.chain") # check the full path

HepG2.hg38 = unlist(liftOver(HepG2, ch))

HepG2.hg19[1]
HepG2.hg38[1]

start(HepG2.hg19[1]) - start(HepG2.hg38[1])
