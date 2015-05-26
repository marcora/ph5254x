library(rafalib)
library(SpikeIn)

data(SpikeIn95)

spms = pm(SpikeIn95)
spd = pData(SpikeIn95)

mypar(1,1)

# box plot of probe intensities
boxplot(log2(spms), xlab="technical replicates", ylab="log2(intensity)", range=0)

# smooth histogram of probe intensities
shist(log2(spms[,1]), unit=0.1, type="n", xlab="log2(intensity)", main="Ten technical replicates")
for (i in 1:10) {
  shist(log2(spms[,i]), unit=0.1, col=i, add=TRUE, lwd=2, lty=i)
}

# M-A plot for arrays 9 and 10
i=10
j=9

M = log2(spms[,i]) - log2(spms[,j])
A = (log2(spms[,i]) + log2(spms[,j]))/2

snames = colnames(spd)[which(spd[i, ]/spd[j, ] == 2)] # probes for spiked-in RNAs with fold-changes of 2
sidx = which(probeNames(SpikeIn95) %in% snames)

splot(A, M, ylim=c(-1.5,1.5))
points(A[sidx], M[sidx], ylim=c(-4,4), pch=21, bg=1)

###############################################################################

# reduce number of data points for faster calculation of LOESS fit

o = order(A) # to sample along the A axis
m = M[o]
a = A[o]

# for a sample of equally spaced data points along the A axis
ind = round(seq(1, length(a), length.out = 5000))

a = a[ind]
m = m[ind]

fit = loess(m ~ a)

mypar(1,1)
splot(A, M, ylim=c(-1.5,1.5))
points(A[sidx], M[sidx], ylim=c(-4,4), pch=21, bg=1)
lines(a, fit$fitted, lwd=2, col=2)

# The default LOESS method in R fits parabulas instead of lines and results in the
# undesired curving up at the right, driven by a few data points with high values of A
# We can avoid the default behavior and force the LOESS method to fit lines by specifying degree=1
fit = loess(m ~ a, degree=1)

lines(a, fit$fitted, lwd=2, col=3)

# The default LOESS method in R fits using 3/4 of the data points.
# We can make the line more "flexible" by specifying a different value for span, e.g., 0.5
fit = loess(m ~ a, degree=1, span=.5)

lines(a, fit$fitted, lwd=2, col=4)

# now plot the LOESS normalized data, i.e., the data with the non-linear effect estimated by LOESS removed
bias = predict(fit, newdata=data.frame(a=A))
nM = M - bias # normalized M

splot(A, nM, ylim=c(-1.5,1.5))
points(A[sidx], nM[sidx], ylim=c(-4,4), pch=21, bg=1)
abline(h=0, lwd=2, col=2)

###############################################################################

# quantile normalization

library(preprocessCore)

nspms = normalize.quantiles(spms)

# M-A plot for arrays 9 and 10 before quantile normalization
i=10
j=9

M = log2(spms[,i]) - log2(spms[,j])
A = (log2(spms[,i]) + log2(spms[,j]))/2

splot(A, M, ylim=c(-1.5,1.5))
points(A[sidx], M[sidx], ylim=c(-4,4), pch=21, bg=1)

# M-A plot for arrays 9 and 10 after quantile normalization
i=10
j=9

M = log2(nspms[,i]) - log2(nspms[,j])
A = (log2(nspms[,i]) + log2(nspms[,j]))/2

splot(A, M, ylim=c(-1.5,1.5))
points(A[sidx], M[sidx], ylim=c(-4,4), pch=21, bg=1)


