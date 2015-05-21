library(affy)

basedir = file.path(getwd(), "rawdata/celfiles")
sampleinfo = read.delim(file.path(basedir, "sampleinfo.txt"), check.names=FALSE, as.is=TRUE)
rownames(sampleinfo) = sampleinfo$filenames
head(sampleinfo)
filenames = list.celfiles(basedir)
filenames %in% sampleinfo[, 1]

ab = ReadAffy(filenames=file.path(basedir, sampleinfo[, 1]), phenoData=sampleinfo)

class(ab)

dim(pm(ab))
dim(pData(ab))
annotation(ab)

rownames(pData(ab))
colnames(pm(ab))
es = rma(ab)

dim(pm(ab))
dim(exprs(es))

es = justRMA(filenames=file.path(basedir, sampleinfo[, 1]), phenoData=sampleinfo)

###############################################################################

# source("http://bioconductor.org/biocLite.R")
# biocLite("limma")

library(limma)

basedir = file.path(getwd(), "rawdata/agilent")
targets = readTargets(file.path(basedir, "TargetBeta7.txt"))
head(targets)

RG = read.maimages(targets$FileNames, source="genepix", path=basedir)

class(RG)
dim(RG$R)
dim(RG$G)

MA = MA.RG(RG, bc.method="none")

class(MA)
dim(MA$M)
dim(MA$A)

# MA plot of first sample
plot(MA$A[,1], MA$M[,1])

# image plot of second sample
imageplot(MA$M[,2], RG$printer, zlim=c(-3,3))
