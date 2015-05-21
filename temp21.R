library(pasillaBamSubset)
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
txdb = TxDb.Dmelanogaster.UCSC.dm3.ensGene

gs = genes(txdb)
gs = gs[seqnames(gs) == "chr4"]
class(gs)

grl = exonsBy(txdb, by="gene")
grl = grl[names(gs)]
class(grl)

all.equal(names(gs), names(grl))


bf = BamFile(untreated1_chr4())

gs.so = summarizeOverlaps(features=gs,
                           reads=bf,
                           ignore.strand=TRUE)

grl.so = summarizeOverlaps(features=grl,
                           reads=bf,
                           ignore.strand=TRUE)

# adding a pseudocount fixes the log(0) issue
plot(assay(gs.so)+1, assay(grl.so)+1,log="xy")
abline(0,1)
ratio = assay(grl.so)/assay(gs.so)
mean(ratio[assay(gs.so) > 0])


count = assay(grl.so)

# here, we can just use sum() because there is only one column
# otherwise we would use: sweep(count, 2, colSums(count), "/")

fpm = (count/sum(count)) * 1e6
head(fpm, 1)

###############################################################################

ebp = sum(width(reduce(grl)))
summary(ebp)

fpkm = fpm/ebp*1e3
head(fpkm, 1)
