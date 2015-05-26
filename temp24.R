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

