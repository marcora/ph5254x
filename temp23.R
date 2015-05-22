# library(devtools)
# install_github("genomicsclass/maPooling")

library(maPooling)
data(maPooling)

dat = c(rep.int(0, 617), rep.int(1, 306), rep.int(2, 65), rep.int(3, 12))
hist(dat)

model.pois = rpois(length(dat), 0.5)
hist(model.pois)

model.binom = rbinom(length(dat), size = 5e6, p = 0.5)
hist(model.binom)

model.norm = rnorm(length(dat))
hist(model.norm)

u = exprs(maPooling)[,1]
v = exprs(maPooling)[,2]
x = exprs(maPooling)[,3]
y = exprs(maPooling)[,4]

cor(u, v)
cor((u-v), (x-y))
cor(log(u/v), log(x/y))

###############################################################################

# source("http://bioconductor.org/biocLite.R")
# biocLite("SpikeIn")
# biocLite("SpikeInSubset")
# biocLite("hgu133atagcdf")

library(rafalib)
library(affy)
library(SpikeIn)
library(SpikeInSubset)
library(hgu133atagcdf)

# data(package="SpikeInSubset")

data(SpikeIn133)
head(pData(SpikeIn133))

colnames(pData(SpikeIn133))

pData(SpikeIn133)[ ,"203508_at"]

###############################################################################

library(SpikeInSubset)
data(mas133)
e=exprs(mas133)##get expression
A=(log2(e[,4])+log2(e[,1]))/2
M=log2(e[,4]/e[,1])

##find the genes that were spiked in to have 
##fold changes of 2
siNames=colnames(pData(mas133))
siNames=siNames[pData(mas133)[4,]/pData(mas133)[1,]==2]
spikeinIndex=match(siNames,rownames(e))

mypar2(1,1)
splot(A,M,ylim=c(-4,4),cex=0.5)
abline(h=c(-1,1),col=1,lwd=2,lty=2)
points(A[spikeinIndex],M[spikeinIndex],bg=2,pch=21)

###############################################################################

data(SpikeIn133)
pd=pData(SpikeIn133)[1:14,] ##pick the first 14, rest are reps
pns=probeNames(SpikeIn133)
pms=pm(SpikeIn133)[,1:14] ##pick the first 14, rest are reps
ind=which(pns==colnames(pd)[1]) ##probes in gene 1
concentration=pd[,1]
concentration[concentration==0]= 1/16
mypar2(1,1)
matplot(log2(concentration),t(log2(pms[ind,])),xlab="log (base 2) concentration",ylab="log (base 2) intensity", type="l")

###############################################################################

pd = pData(SpikeIn133) # sample meta data
pms = pm(SpikeIn133) # perfect-match probe intensities
pns = probeNames(SpikeIn133) # probe names

j = which(colnames(pd) == "203508_at")
concentration=pd[,j]
col.idx = which(concentration==0)

row.idx = which(pns == "203508_at")

min(pms[row.idx,col.idx])
max(pms[row.idx,col.idx])

###############################################################################

mms = mm(SpikeIn133) # mis-match probe intensities

y = as.vector(pms[row.idx,col.idx])
x = as.vector(mms[row.idx,col.idx])
plot(x,y)
plot(log2(x),log2(y))
cor(log2(x),log2(y))

hist(log2(x)-log2(y))

###############################################################################

bg1 = bg.correct.mas(SpikeIn133)
bg2 = bg.correct.rma(SpikeIn133)

pd= pData(SpikeIn133)
pns=probeNames(SpikeIn133)
pms1=pm(bg1) 
pms2=pm(bg2)

ind=which(pns==colnames(pd)[1]) ##probes in gene 1
concentration=pd[,1]
concentration[concentration==0]= 1/16

mypar2(1,2)
matplot(log2(concentration),t(log2(pms1[ind,])),xlab="log (base 2) concentration",ylab="log (base 2) instensity",ylim=c(0,13))
matplot(log2(concentration),t(log2(pms2[ind,])),xlab="log (base 2) concentration",ylab="log (base 2) instensity",ylim=c(0,13))

pData(SpikeIn133)[c(1,15,29),]
pData(SpikeIn133)[c(2,16,30),]

ind = c(1,15,29)
pm1 = log2( pm(bg1)[,ind])
pm2 = log2( pm(bg2)[,ind])

SD1 = rowSds(pm1)
A1 = rowMeans(pm1)
SD2 = rowSds(pm2)
A2 = rowMeans(pm2)
mypar2(2,1)
splot(A1,SD1,ylim=c(0,3),cex=.25)
splot(A2,SD2,ylim=c(0,3),cex=.25)