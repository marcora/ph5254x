###############################################################################
## ExpressionSet (for microarray data)                                       ##
###############################################################################

# source("http://bioconductor.org/biocLite.R")
# biocLite("GEOquery")

library(Biobase)

library(GEOquery)

geoq = getGEO("GSE9514")

e = geoq[[1]]

e

dim(e)
nrow(e)
ncol(e)

# expression data as a matrix (probesets intensities x samples)
exprs(e)[1:3, 1:3]

# meta data/information about the columns/samples
pData(e)

dim(pData(e))

names(pData(e))

pData(e)$characteristics_ch1 # conditions

table(pData(e)$characteristics_ch1) # to inspect the number of replicates for each condition

# meta data/information about the rows/probesets
fData(e)

dim(fData(e))

names(fData(e))

head(rownames(fData(e)))

# meta data/information about the experiment and other annotation
experimentData(e)
annotation(e)

###############################################################################
## SummarizedExperiment (for ngs data)                                       ##
###############################################################################

# source("http://bioconductor.org/biocLite.R")
# biocLite("parathyroidSE")

library(Biobase)

library(parathyroidSE)

data(parathyroidGenesSE)

se = parathyroidGenesSE

se

dim(se)
nrow(se)
ncol(se)

# expression data as a matrix (counts x samples)
assay(e)[1:3, 1:3]

# meta data/information about the columns/samples
colData(se)

names(colData(se))

colData(se)$treatment # conditions

table(colData(se)$treatment) # to inspect the number of replicates for each condition

# meta data/information about the rows/genomic ranges -> a GRangesList
rowData(se)

length(rowData(se)) # total number of genes

rowData(se)[[1]] # GRanges object containing information about the exons of the first gene
length(rowData(se)[[1]]) # total number of exons

head(rownames(se)) # gene identifiers
length(rownames(se)) # total number of genes

# meta data/information about the experiment and other annotation
metadata(rowData(se))

exptData(se)
exptData(se)$MIAME

