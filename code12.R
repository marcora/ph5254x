library(Rsamtools)
library(GenomicRanges)
library(GenomicAlignments)

library(pasillaBamSubset)
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)

txdb = TxDb.Dmelanogaster.UCSC.dm3.ensGene
grl = exonsBy(txdb, by="gene")
grl[100] # => GRangesList
grl[[100]] # => GRanges object for gene 100
grl[[100]][1] # => GRanges object for exon 1 of gene 100

fl1 = untreated1_chr4()
fl2 = untreated3_chr4()

fls = BamFileList(c(fl1, fl2)) # specify yieldSize (e.g., 10e6) when reading real BAM files, depending on available RAM!!!
names(fls) = c("first","second")

so1 = summarizeOverlaps(features=grl,
                        reads=fls,
                        ignore.strand=TRUE)

so1

head(assay(so1))
colSums(assay(so1))

rowData(so1)
colData(so1)
colData(so1)$sample = c("one","two")
colData(so1)

metadata(rowData(so1)) 

x = assay(so1)[,1]
hist(x[x > 0], col="grey")
hist(x[x > 0 & x < 10000], col="grey")
plot(assay(so1) + 1, log="xy")


fls = BamFileList(fl2)
so2 = summarizeOverlaps(features=grl,
                        reads=fls,
                        ignore.strand=TRUE,
                        singleEnd=FALSE,
                        fragments=TRUE)
head(assay(so2))
colSums(assay(so2))

rowData(so2)
colData(so2)
colData(so2)$sample = c("one","two")
colData(so2)

metadata(rowData(so2)) 

x = assay(so2)[,1]
hist(x[x > 0], col="grey")
hist(x[x > 0 & x < 10000], col="grey")
plot(assay(so1) + 1, log="xy")


plot(assay(so1)[,2], assay(so2)[,1], xlim=c(0,5000), ylim=c(0,5000),
     xlab="single end counting", ylab="paired end counting")
abline(0,1)
abline(0,.5)
