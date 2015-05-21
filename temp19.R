library(hgu95acdf)

library(affy)

basedir = file.path(getwd(), "rawdata/celfiles")

sampleinfo = read.delim(file.path(basedir, "sampleinfo.txt"), check.names=FALSE, as.is=TRUE)
rownames(sampleinfo) = sampleinfo$filenames
head(sampleinfo)

sampleinfo["1521a99hpp_av06.CEL.gz","36311_at"]

filenames = list.celfiles(basedir)
filenames %in% sampleinfo[, 1]

ab = ReadAffy(filenames=file.path(basedir, sampleinfo[, 1]), phenoData=sampleinfo)

sum(probeNames(ab) == '36311_at')

length(featureNames(ab))
length(probeNames(ab))

mat = pm(ab)[probeNames(ab) == '36085_at', ]
mat

conc = pData(ab)[c('1532a99hpp_av04.CEL.gz', '1532b99hpp_av04.CEL.gz'), '36085_at']
conc

par(mfrow=c(2,1))
matplot(mat[, c('1532a99hpp_av04.CEL.gz', '1532b99hpp_av04.CEL.gz')], type="l", log="y")
plot(log2(mat[, '1532b99hpp_av04.CEL.gz']/mat[, '1532a99hpp_av04.CEL.gz']), type="l")
abline(h=log2(conc[2]/conc[1]))
par(mfrow=c(1,1))

###############################################################################

es = rma(ab)

g = as.factor(c(1,1,1,2,2,2))

library(genefilter)

tt = rowttests(es, g)
tt['36085_at', 'p.value']

###############################################################################

sig = colnames(pData(ab))[-1]
sig

tt$spiked = rownames(tt) %in% sig

boxplot(-dm ~ spiked, data = tt, xlab = "spiked")
boxplot(-dm ~ spiked, data = tt, xlab = "spiked", ylim=c(-1,1))

###############################################################################

library(limma)

basedir = file.path(getwd(), "rawdata/agilent")
targets = readTargets(file.path(basedir, "TargetBeta7.txt"))
head(targets)

RG = read.maimages(targets$FileNames, source="genepix", path=basedir)

head(RG$genes)

i = which(RG$genes$ID=="H200015482")
j = which(rownames(RG$targets)=="6Hs.168")
log2(RG$R[i,j]/RG$G[i,j])

# or

MA = MA.RG(RG,bc.method="none")
MA$M[i,j]

###############################################################################

