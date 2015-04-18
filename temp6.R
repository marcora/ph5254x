# library(devtools)
# install_github("genomicsclass/ph525x")
library(GenomicRanges)
library(ph525x) # for plotRanges

gr = GRanges("chrZ", IRanges(start=c(5,10), end=c(35,45)), "+", seqlengths=c(chrZ=100L))
gr

shift(gr, 10)
shift(gr, 80)
trim(shift(gr, 80))
mcols(gr)
mcols(gr)$value = c(-1,4)
mcols(gr)$value = NULL

gr2 = GRanges("chrZ", IRanges(11:13, 51:53))
gr2

# use grangeslist to, e.g., group exons by gene/transcript
grl = GRangesList(gr, gr2)
grl

length(grl) # a list with two elements/granges

gr1 = GRanges("chrZ", IRanges(c(1,11,21,31,41), width=5))
gr2 = GRanges("chrZ", IRanges(c(19,33), c(38,35)))
gr1
gr2

fo = findOverlaps(gr1, gr2)
fo
queryHits(fo)
subjectHits(fo)
gr1 %over% gr2

# see also Rle and Views(Rle, ...)


x = GRanges("chr1", IRanges(c(1,101),c(50,150)), strand=c("+","-"))
x

plotGRanges = function(x, ...) plotRanges(ranges(x), ...)

par(mfrow=c(2,1))
plotGRanges(x)
plotGRanges(resize(x,1))



# simulate two diff gene transcript isoforms
x = GRanges("chr1", IRanges(c(101,201,401,501),c(150,250,450,550)), strand="+")
y = GRanges("chr1", IRanges(c(101,221,301,401,541),c(150,250,350,470,550)), strand="+")
x
y

par(mfrow=c(2,1))
plotGRanges(x)
plotGRanges(y)

grl = GRangesList(x,y)
grl

xy = c(x,y)
xy

par(mfrow=c(5,1))
plotGRanges(xy, main="original", xlim=c(0,600))
plotGRanges(disjoin(xy), main="disjoin", xlim=c(0,600))
plotGRanges(reduce(xy), main="reduce", xlim=c(0,600))
plotGRanges(intersect(x,y), main="intersect", xlim=c(0,600))
plotGRanges(setdiff(reduce(xy), intersect(x,y)), main="!intersect", xlim=c(0,600))

sum(width(xy))
sum(width(disjoin(xy)))
sum(width(reduce(xy)))
sum(width(intersect(x,y)))
sum(width(setdiff(reduce(xy), intersect(x,y))))

disjoined = disjoin(c(x,y))
in.both = disjoined %over% x & disjoined %over% y
sum(width(disjoined[in.both]))

not.in.both = !(disjoined %over% x & disjoined %over% y)
sum(width(disjoined[not.in.both]))

z = GRanges("chr1", range(ranges(x)), "-")
z

sum(x %over% z)
