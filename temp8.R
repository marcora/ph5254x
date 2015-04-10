library(Homo.sapiens)

ghs = genes(Homo.sapiens)

genome(ghs)
length(ghs)
sort(table(seqnames(ghs)))

hist(width(ghs))
hist(log10(width(ghs)))
median(width(ghs))
