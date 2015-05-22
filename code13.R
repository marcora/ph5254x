# library(devtools)
# install_github("genomicsclass/SpikeInEDA")

library(SpikeInEDA)
data(SpikeInEDA)

hist(int[,1])
hist(log2(int[,1]))

# overlayed density plots, microarray #4 is problematic!
for (i in 1:ncol(int)) {
  if (i==1) {
    plot(density(log2(int[,i])), col=(i==4)+1)
  } else {
    lines(density(log2(int[,i])), col=(i==4)+1)
  }
}

# everything is highly correlated! Correlation is not good at identifying a problematic microarray
signif(cor(int), 2)

# scatter plots (n.b., use a random subset of all the points)

splot <- function(x, y, ...) {
  ind = sample(length(x), 10000)
  x = x[ind]
  y = y[ind]
  plot(x, y, ...)
}

par(mfrow=c(2,1))
splot(log2(int[,1]), log2(int[,2]))
splot(log2(int[,1]), log2(int[,4]))

# ma plots (n.b., use a random subset of all the points)

maplot <- function(x, y, ...) {
  splot((x+y)/2, y-x, ...)
}

par(mfrow=c(2,1))
maplot(log2(int[,1]), log2(int[,2]))
maplot(log2(int[,1]), log2(int[,4]))

###############################################################################

library(matrixStats)
i = 4 # or 1
r = log2(int[,i]) - rowMedians(log2(int))

MAX = 1 # to avoid outliers taking over colors of image
r[r>MAX] = MAX
r[r<-MAX] = -MAX

mat = matrix(NA, max(locations[,1]), max(locations[,2]+1)/2)
for (j in 1:nrow(locations)) {
  mat[locations[j,1], (locations[j,2]+1)/2] = r[j]
}

par(mfrow=c(1,1))
image(mat, col=RColorBrewer::brewer.pal(11, 'RdBu'))
