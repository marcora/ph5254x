# library(devtools)
# install_github("genomicsclass/ERBS")
library(ERBS)

data(HepG2)
data(GM12878)

HepG2 # HepG2 cell line
GM12878 # GM12878 cell line

# overlaps
fo = findOverlaps(HepG2, GM12878)
erbs = HepG2[queryHits(fo)]
erbs = granges(erbs)

# intersection
erbs2 = intersect(HepG2, GM12878)
erbs2 = granges(erbs2)

# sort
erbs = sort(erbs)
erbs2 = sort(erbs2)

# confirm same chr and strand
all(seqnames(erbs) == seqnames(erbs2))
all(strand(erbs) == strand(erbs2))

# proportion of identical ranges based on start and end
mean(start(erbs) == start(erbs2) & end(erbs2) == end(erbs3))

# the intersection should be smaller
all(width(erbs2) <= width(erbs))

###############################################################################
library(Homo.sapiens)

ghs = genes(Homo.sapiens)

tss = resize(ghs, 1)

start(tss["100113402"])

fo = findOverlaps(HepG2, GM12878)
erbs = HepG2[queryHits(fo)]

index = nearest(erbs[4], tss)

geneid = names(tss[index])
geneid

keytypes(Homo.sapiens)
columns(Homo.sapiens)

select(Homo.sapiens, keytype="GENEID", keys=geneid, columns = c("SYMBOL"))
