# library(devtools)
# install_github("genomicsclass/ph525x")
library(ph525x)
f1 = dir(system.file("extdata", package="ERBS"), full=TRUE)[1]
f1

readLines(f1, 4)
cat(readLines(f1, 4), sep="\n")


library(tracklayer)
imp = import(f1, format="bedGraph")
class(imp)
imp

genome(imp) = "hg19"
metadata(imp) = list(celltype = "Immortalized B cell")
imp

export(imp, "demoex.bed")
cat(readLines("demoex.bed", 4), sep="\n")
