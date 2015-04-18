# library(devtools)
# install_github("genomicsclass/ERBS")
library(ERBS)

data(HepG2)
data(GM12878)

HepG2 # HepG2 cell line
GM12878 # GM12878 cell line

res = findOverlaps(HepG2, GM12878)
res

index = queryHits(res)

erbs = HepG2[index]
erbs

erbs = granges(erbs)
erbs

# erbs contains the ER binding sites in HepG2 that overlap with those in GM12878

# which genes contain ER binding sites in their promoter?

library(Homo.sapiens)

ghs = genes(Homo.sapiens)
ghs

index = precede(erbs, ghs)
index

ghs[index[1:3]]
erbs[1:3]

###############################################################################

# which genes contain ER binding sites near their tss?

tsss = resize(ghs, 1)
tsss

ds = distanceToNearest(erbs, tsss)
ds

dists = mcols(ds)$distance
hist(log10(dists))

tsss_nearest_to_erbs = subjectHits(ds)[dists < 1000]

keytypes(Homo.sapiens)
columns(Homo.sapiens)

keys = as.character(mcols(ghs[index])$GENEID)

select(Homo.sapiens, keytype="GENEID", keys=keys, columns = c("SYMBOL"))
