library(Biobase)
load("bottomly_eset.RData")

bottomly.eset

head(exprs(bottomly.eset))

library(rafalib)

mypar(1,1)

# Examine the histogram of sample-sample correlations of log counts:

hist(cor(log(exprs(bottomly.eset) + 1)))


# There are a number of reasons for such high correlations:
  
# Biological correlation: Genes are expressed similarly across organisms
# Systematic bias: Technical bias affects the same genes similarly (for example, if the bias is due to the sequence of the gene, then this could be common across experiments)
# The presence of many zeros, which anchor the cloud of log counts at (0,0)

# Remove some of the rows with many zeros, and plot the histogram of
# correlations again:

mat = exprs(bottomly.eset)
rs = rowSums(mat)
hist(cor(log(mat[rs>10,] + 1)))

e = bottomly.eset[ rowSums(exprs(bottomly.eset)) > 2, ]

dim(e)

# Calculate the GC content for these mouse genes

library(GenomicFeatures)

source("http://bioconductor.org/biocLite.R")
biocLite("TxDb.Mmusculus.UCSC.mm9.knownGene")
biocLite("BSgenome.Mmusculus.UCSC.mm9")
biocLite("org.Mm.eg.db")

library("TxDb.Mmusculus.UCSC.mm9.knownGene")
library("BSgenome.Mmusculus.UCSC.mm9")
library("org.Mm.eg.db")

res = select(org.Mm.eg.db, keys=rownames(e), keytype="ENSEMBL", columns="ENTREZID")

fData(e)$ENTREZ = res$ENTREZID[ match(rownames(e), res$ENSEMBL) ]

sum(is.na(fData(e)$ENTREZ))
e = e[!is.na(fData(e)$ENTREZ),]
