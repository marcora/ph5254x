library(TxDb.Hsapiens.UCSC.hg19.knownGene)

txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
txdb

exbg = exonsBy(txdb, by="gene")
exbg
