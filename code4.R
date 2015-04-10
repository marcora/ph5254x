# library(devtools)
# install_github("genomicsclass/ERBS")
library(ERBS)

data(HepG2)
data(GM12878)

HepG2 # HepG2 cell line
GM12878 # GM12878 cell line

library(BSgenome.Hsapiens.UCSC.hg19)

Hsapiens

chr17 = Hsapiens$chr17
chr17

hepseq = getSeq(Hsapiens, HepG2)
hepseq

length(hepseq) == length(HepG2)
all(width(hepseq) == width(HepG2))

rhepseq = getSeq(Hsapiens, shift(HepG2, 1000))

mot = "TCAAGGTCA"

sum(vcountPattern(mot, hepseq)) + sum(vcountPattern(mot, reverseComplement(hepseq)))
sum(vcountPattern(mot, rhepseq)) + sum(vcountPattern(mot, reverseComplement(rhepseq)))
