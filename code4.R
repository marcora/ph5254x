library(Biostrings)
library(rtracklayer)

# library(devtools)
# install_github("genomicsclass/ERBS")
library(ERBS)

available.genomes()

library(BSgenome.Hsapiens.UCSC.hg19)

Hsapiens


grep("Drerio", available.genomes(), value=TRUE) # exclude masked


library(BSgenome.Hsapiens.UCSC.hg19.masked)
c17m = BSgenome.Hsapiens.UCSC.hg19.masked$chr17
class(c17m)

c22m = BSgenome.Hsapiens.UCSC.hg19.masked$chr22
c22m

ch = import.chain("hg38ToHg19.over.chain")
ch


library(TxDb.Hsapiens.UCSC.hg38.knownGene)

g38 = genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
