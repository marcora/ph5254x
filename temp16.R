library(Homo.sapiens)
g = genes(Homo.sapiens)
library(ERBS)
data(HepG2)

kp = g[resize(g,1) %over% HepG2]


# source("http://bioconductor.org/biocLite.R")
# biocLite("ReportingTools")

nn = names(kp)
m = select(Homo.sapiens, keys=nn, keytype="ENTREZID",
           columns=c("SYMBOL", "GENENAME", "TERM", "GO"))
library(ReportingTools)
hrep = HTMLReport(shortName="erhep.html")
publish(m, hrep)
finish(hrep)
