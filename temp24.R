# source("http://bioconductor.org/biocLite.R")
# biocLite("SpikeIn")

library(genefilter)
library(SpikeIn)
data(SpikeIn95)

snames = colnames(pData(SpikeIn95))
sidx = which(probeNames(SpikeIn95) %in% snames)

# the spiked-in probes change from array to array with the exception of the last eight arrays
pData(SpikeIn95)[52:59,]

# for the purpose of this assessment, we will remove the arrays in which we have repeated spike-in concentrations
pm = pm(SpikeIn95)[,c(1:51,52,56)]

# we can see the need for normalization by looking at boxplots  
boxplot(log2(pm),range=0)

library(matrixStats)

sds = rowSds(log2(pm))

(median(sds[-sidx])) # first summary (non spiked-in genes)
(median(sds[sidx])) # second summary (spiked-in genes)

boxplot(sds[-sidx],sds[sidx],range=0)

# Note that the first summary relates to specificity: we want measurements
# across replicate arrays to be similar and thus the variability low. However,
# because we designed the experiment for concentrations in the spike-ins to vary
# across samples, we want the concentration of spike-in genes to change so that
# variability should be higher. A successful normalization approach will improve
# specificity (lower the SD of non-spiked-in genes) without affecting
# sensitivity (leave SD of spiked-in genes the same).

# Now use the normalize.quantiles in the preprocessCore package to normalize all
# the probes together (make sure to normalize in the original scale)

library(preprocessCore)

npm = normalize.quantiles(pm)

nsds = rowSds(log2(npm))

(median(nsds[-sidx])) # first summary (non spiked-in genes)
(median(nsds[sidx])) # second summary (spiked-in genes)

boxplot(nsds[-sidx],nsds[sidx],range=0)

# So far, we used the variance of the non-spiked in probes across replicates as
# a measure of specificity and the variance across spiked-in probes as measure a
# sensitivity. A more sophisticated measure of sensitivity involves computing 
# the slope obtained from fitting a line to the the observed intensities versus 
# intended concentrations plot (in the log scale). This slope should be 1 since 
# we expect the observed intensities to double when we add twice as many 
# transcripts.

###############################################################################

# Normalization techniques such as quantile normalization are not always
# appropriate. An example dataset comes from the dataset described in the paper
# Loven et al. (2012)

# library(devtools)
# install_github("stephaniehicks/mycAffyData")

# source("http://bioconductor.org/biocLite.R")
# biocLite("primeviewcdf")
# biocLite("SQN")
# ?SQN

library(mycAffyData)
data(mycData)
erccIndex=grep("ERCC",probeNames(mycData))

# The ERCC stands for the external RNA control consortium, which is trying to
# define standards for high throughput technologies.

# Let's consider three approaches: 1) no normalization, 2) quantile
# normalization, 3) sqn normalization.

library(preprocessCore)
library(SQN)

pms = pm(mycData)
pms = list( log2(pms ),
            log2( normalize.quantiles( pms )),
            SQN(log2(pms),ctrl.id=erccIndex))
names(pms)=c("none","qn","sqn")

# This is an experiment in which we expect most genes to be expressed higher in
# samples 3,4 compared to 1,2. The control genes should not change. To explore
# which normalization technique worked best, compute the difference between
# samples 1 and 3.

M = sapply(pms, function(y){
  y[,3]-y[,1]
})

head(M)

library(rafalib)

mypar(1,2)
boxplot(M[erccIndex,], use.cols=TRUE, main="spiked-in controls", ylim=c(-6,6), range=0)
abline(h=0,lty=2)
boxplot(M[-erccIndex,], use.cols=TRUE, ylim=c(-6,6), range=0)
abline(h=0,lty=2)

