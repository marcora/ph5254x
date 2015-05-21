library(pasillaBamSubset)
library(Rsamtools)

filename = untreated1_chr4()

bf = BamFile(filename)

seqinfo(bf)
quickBamFlagSummary(bf)

gr = GRanges("chr4", IRanges(440000, 470000))

countBam(bf, param=ScanBamParam(which=gr))
# or
reads.pos = scanBam(bf, param=ScanBamParam(what="pos", which=gr))
length(reads.pos[[1]]$pos)

###############################################################################

library(Biostrings)

reads.seq = scanBam(bf, param=ScanBamParam(what="seq", which=gr))

mean(letterFrequency(reads[[1]]$seq, "GC", as.prob=TRUE))

###############################################################################

ga = readGAlignments(BamFile(filename))

hist(start(ga), breaks=100)

library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)

gs = genes(TxDb.Dmelanogaster.UCSC.dm3.ensGene)

gs2 = gs[gs %over% GRanges("chr4", IRanges(200000, 300000))]

head(gs2)

countOverlaps(gs2['FBgn0039890'], ga)
