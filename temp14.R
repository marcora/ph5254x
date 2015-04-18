library(rtracklayer)

data(targets)

class(targets)

head(targets)

library(GenomicRanges)

mtar = with(targets,
            GRanges(chrom, IRanges(start,end), strand=strand,
                    targets=target, mirname=name))
mtar

cat(export(mtar[1:5], format="bed"), sep="\n")
cat(export(mtar[1:5], format="gff3"), sep="\n")
