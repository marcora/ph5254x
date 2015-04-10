# library(devtools)
# install_github("genomicsclass/ph525x")

library(IRanges)

ir = IRanges(101,200)

ir * 2 # zoom in or out (-2) by multipling/dividing the width by the given factor

narrow(ir, start=20) # add/remove from one side or another

ir + 25 # add to both sides

ir - 25 # remove from both sides

ir = IRanges(start=c(1,11,21), end=c(3,15,27))
ir

width(ir)
sum(width(ir))

ir = IRanges(start=c(101,106,201,211,221,301,306,311,351,361,401,411,501), end=c(150,160,210,270,225,310,310,330,390,380,415,470,510))
ir

library(ph525x)

par(mfcol=c(5,1))
plotRanges(ir, xlim=c(0,600))

range_ir = range(ir)
plotRanges(range_ir, main="range", xlim=c(0,600))

reduce_ir = reduce(ir)
plotRanges(reduce_ir, main="reduce", xlim=c(0,600))

gaps_ir = gaps(ir)
plotRanges(gaps_ir, main="gaps", xlim=c(0,600))

disjoin_ir = disjoin(ir)
plotRanges(disjoin_ir, main="disjoin", xlim=c(0,600))


# What is the total width from 101 to 510 which is not covered by ranges in ir?
sum(width(gaps(ir)))

# How many disjoint ranges are contained within ir?
length(disjoin(ir))


par(mfrow=c(4,1))

plotRanges(ir, xlim=c(0,600))
plotRanges(resize(ir,1), xlim=c(0,600)) # fix="start" is the default
plotRanges(resize(ir,1, fix="center"), xlim=c(0,600))
plotRanges(resize(ir,1, fix="end"), xlim=c(0,600))

###############################################################################

library(GenomicRanges)

gr <- GRanges(seqnames="chr1", strand="+", ranges=IRanges(start=c(1,3,5), end=c(3,5,7)))
gr

# accessor functions
length(gr)
names(gr)
seqnames(gr)
ranges(gr)
strand(gr)
mcols(gr)
seqlengths(gr)
seqinfo(gr)
genome(gr)
isCircular(gr)
score(gr)
