###############################################################################
## ExpressionSet (for microarray data)                                       ##
###############################################################################

# source("http://bioconductor.org/biocLite.R")
# biocLite("GEOquery")

library(Biobase)

library(GEOquery)

geoq = getGEO("GSE9514")

es = geoq[[1]]

es

dim(es)
nrow(es)
ncol(es)

# expression data as a matrix (probesets intensities x samples)
exprs(es)[1:3, 1:3]

# metadata/information about the columns/samples
pData(es)

dim(pData(es))

names(pData(es))

pData(es)$characteristics_ch1 # conditions

table(pData(es)$characteristics_ch1) # to inspect the number of replicates for each condition

# metadata/information about the rows/probesets
head(fData(es))

dim(fData(es))

names(fData(es))

head(rownames(fData(es)))

# meta data/information about the experiment and other annotation
experimentData(es)
annotation(es)

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
assay(se)[1:3, 1:3]

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
