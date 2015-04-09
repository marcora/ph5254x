# source("http://bioconductor.org/biocLite.R")
# or
# library(BiocInstaller)
# biocLite("Homo.sapiens")

library(tidyr)
library(dplyr)
library(magrittr)

library(Homo.sapiens)
class(Homo.sapiens)

keytypes(Homo.sapiens)
columns(Homo.sapiens)

head(keys(Homo.sapiens, keytype="ENTREZID"))
length(unique(keys(Homo.sapiens, keytype="ENTREZID")))

head(keys(Homo.sapiens, keytype="ENSEMBL"))
length(unique(keys(Homo.sapiens, keytype="ENSEMBL")))

select(Homo.sapiens, key="9575", keytype="ENTREZID", columns=c("SYMBOL", "GENENAME", "ENSEMBL", "ENTREZID"))

tab = select(Homo.sapiens, key="circadian rhythm", keytype="TERM", columns=c("SYMBOL", "GENENAME", "ENSEMBL", "ENTREZID"))

length(unique(tab$ENTREZID))
