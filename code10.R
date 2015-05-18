library(affy)

basedir = file.path(getwd(), "rawdata/celfiles")
tbl = read.delim(file.path(basedir, "sampleinfo.txt"))
head(tbl)
filenames = list.celfiles(basedir)
filenames %in% tbl[, 1]

ab = ReadAffy(filenames=file.path(basedir, tbl[, 1]), phenoData=tbl)

