# source("http://bioconductor.org/biocLite.R")
# biocLite("pasillaBamSubset")
# biocLite("Rsamtools")

library(pasillaBamSubset)
library(Rsamtools)

(filename = untreated1_chr4())

(bf = BamFile(filename))

seqinfo(bf)

(sl = seqlengths(bf))

quickBamFlagSummary(bf)

(gr = GRanges("chr4", IRanges(1, sl["chr4"])))

countBam(bf, param=ScanBamParam(which = gr))

(reads = scanBam(BamFile(filename, yieldSize=5)))

names(reads[[1]])

gr = GRanges("chr4",IRanges(500000, 700000))
reads = scanBam(bf, param=ScanBamParam(what=c("pos","strand"), which=gr))
hist(reads[[1]]$pos)

readsByStrand = split(reads[[1]]$pos, reads[[1]]$strand)
myHist = function(x) table(cut(x, 50:70 * 10000))
tab = sapply(readsByStrand, myHist)
barplot(t(tab))

###############################################################################

# source("http://bioconductor.org/biocLite.R")
# biocLite("GenomicAlignments")

library(GenomicAlignments)

(ga <- readGAlignments(bf))

gr = GRanges("chr4", IRanges(700000, 800000))
(fo = findOverlaps(ga, gr)) # which reads over this range

countOverlaps(gr, ga) # count overlaps of range with the reads

table(ga %over% gr) # logical vector of read overlaps with the range
