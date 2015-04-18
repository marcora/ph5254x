library(rtracklayer)

ch = import.chain("hg38ToHg19.over.chain")
ch

library(TxDb.Hsapiens.UCSC.hg38.knownGene)

g38 = genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
g38

g19 = liftOver(g38, ch) # return a GRangesList
g19

u19 = unlist(g19)
u19
