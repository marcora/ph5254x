# library(devtools)
# install_github("genomicsclass/ERBS")

library(ERBS)

data(HepG2)
data(GM12878)

HepG2
GM12878

length(HepG2)

# genomic coordinates components
seqnames(HepG2) # chr
ranges(HepG2) # start:end
strand(HepG2) # strand
# metadata columns
mcols(HepG2)
# or
values(HepG2)

seqlengths(HepG2)
names(HepG2)

width(HepG2)
start(HepG2)
end(HepG2)

chrs = seqnames(HepG2)
as.character(chrs)
table(chrs)[1:24]

median(HepG2$signalValue)
seqnames(HepG2[which.max(HepG2$signalValue)])

length(HepG2[seqnames(HepG2) == "chr16"])

hist(width(HepG2))
median(width(HepG2))

###############################################################################

library(IRanges)

ir = IRanges(5,10)
ir

ir = IRanges(5, width=6)
ir

start(ir)
end(ir)
width(ir)

ir = IRanges(start=c(3,5,17), end=c(10,8,20))
ir

length(ir)
start(ir)
end(ir)
width(ir)

ir = IRanges(5,10)
shift(ir, -2)
narrow(ir, start=2)
# see also the flank() function