library(GenomicFeatures)
library(GenomicRanges)
library(IRanges)

# library(devtools)
# install_github("genomicsclass/ERBS")

library(ERBS)

data(HepG2)
data(GM12878)

HepG2
GM12878

start(HepG2[17])

d = distanceToNearest(HepG2[17], GM12878)
i = subjectHits(d)
start(GM12878[i])

d = distanceToNearest(HepG2[17], GM12878)
mcols(d)$distance

ds = distanceToNearest(HepG2, GM12878)

mean(mcols(ds)$distance < 2000)
