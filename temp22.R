if (!file.exists("bottomly_eset.RData")) download.file("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData","bottomly_eset.RData")

load("bottomly_eset.RData")

library("Biobase")

bottomly.eset
dim(bottomly.eset)
pData(bottomly.eset)

p0s = colMeans(exprs(bottomly.eset) == 0)
median(p0s)

###############################################################################

# proportion of genes with zero mapped reads in the three experimental batches
boxplot(split(p0s , pData(bottomly.eset)$experiment.number))


# To further explore batch effects we want to make a multidimensional scaling
# plot. However, when computing distances we prefer to have measures with
# similar variability across samples. Because this is RNAseq data, the variance
# depends on the mean.
library(matrixStats)
A = rowMeans(exprs(bottomly.eset))
SD = rowSds(exprs(bottomly.eset))
plot(A,SD)

# The log transformation is a natural choice here for variance stabilization,
# but we have many 0s so we can't just apply it.
y = log2(exprs( bottomly.eset )+0.5)

library(rafalib)
mypar(2,1)
hist(y[,1],nc=100)
hist(y[y[,1]>0,1],nc=100)
abline(v=3) # = 2^3 = 8 mapped reads

Y = exprs(bottomly.eset)
ind = which(apply(Y >= 8, 1, all))
Y = log2(Y[ind,])

batch = pData(bottomly.eset)$experiment.number
strain = pData(bottomly.eset)$strain

s = svd(Y - rowMeans(Y))
mypar(1,1)
plot(s$d^2/sum(s$d^2), ylab='Variance explained')
d = dist(t(Y))
mds = cmdscale(d)
plot(mds, col = as.fumeric(batch), pch = as.fumeric(strain)+15, xlab = "Dimension 1", ylab = "Dimension 2")
legend("topleft", col=unique(as.fumeric(batch)), legend=unique(batch), pch=1)
legend("bottomleft", pch=unique(as.fumeric(strain)+15), legend=unique(strain))
