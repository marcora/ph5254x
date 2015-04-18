library(Biostrings)

available.genomes()

library(BSgenome.Hsapiens.UCSC.hg19)

Hsapiens


grep("Drerio", available.genomes(), value=TRUE) # exclude masked


library(BSgenome.Hsapiens.UCSC.hg19.masked)
c17m = BSgenome.Hsapiens.UCSC.hg19.masked$chr17
class(c17m)

c22m = BSgenome.Hsapiens.UCSC.hg19.masked$chr22
c22m
