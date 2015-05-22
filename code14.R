if (!file.exists("bottomly_eset.RData")) download.file("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData","bottomly_eset.RData")

load("bottomly_eset.RData")

library("Biobase")

bottomly.eset
dim(bottomly.eset)
pData(bottomly.eset)
head(exprs(bottomly.eset))

# add 0.5 such that the log-transform is -1 for genes with zero reads instead of -Inf

Y = log2(exprs(bottomly.eset)+0.5)

library("rafalib")
mypar(1,1)

for (i in 1:ncol(Y)) {
  shist(Y[,i], unit=0.25, col=i, plotHist=FALSE, add=i!=1)
}

for (i in 1:ncol(Y)) {
  idx = Y[,i] > -1
  shist(Y[idx,i], unit=0.25, col=i, plotHist=FALSE, add=i!=1)
}

mypar(1,2)
idx = rowSums(Y[,1:2]) > 0
plot(Y[idx,1], Y[idx,2], cex=.1)
rm = rowMeans(2^Y[idx,1:2]) # Y is log2(mapped reads)
simulated1 = rpois(length(idx), rm)
simulated2 = rpois(length(idx), rm)
plot(log2(simulated1 + .5), log2(simulated2 + .5), cex=.1)

mypar(1,1)
maplot(Y[idx,1],Y[idx,2])
