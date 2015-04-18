# library(devtools)
# install_github("genomicsclass/GSE5859Subset")

library(GSE5859Subset)
data(GSE5859Subset)
dim(geneExpression)
dim(sampleInfo)
dim(geneAnnotation)

identical(colnames(geneExpression), sampleInfo$filename)
identical(rownames(geneExpression), geneAnnotation$PROBEID)

pd = AnnotatedDataFrame(sampleInfo)

rownames(pd) = colnames(geneExpression)

pData(pd)["GSM136530.CEL.gz","date"]

varLabels(pd)[1]


fd = AnnotatedDataFrame(geneAnnotation)

rownames(fd) = rownames(geneExpression)

pData(fd)["204810_s_at","CHR"]

eset = ExpressionSet(assayData = geneExpression, phenoData = pd, featureData = fd)

annotation(eset) = "hgfocus"

dim(pData(eset))
dim(exprs(eset))
dim(featureData(eset))

ind1 = which(featureData(eset)$CHR=="chrY")
ind2 =  pData(eset)$group==1
femaleY = colMeans(exprs(eset)[ind1, ind2]) 
maleY = colMeans(exprs(eset)[ind1, !ind2]) 
boxplot(maleY,femaleY)
median(maleY)-median(femaleY)


eset[1:10,1:5]

###############################################################################

library(Homo.sapiens)
genes = genes(Homo.sapiens)

library(hgfocus.db) # microarray platform used to obtain the data

# get the ENTREZIDs associated with the PROBEIDs for this array
map = select(hgfocus.db, keys=featureNames(eset), columns="ENTREZID", keytype="PROBEID")

# since we obtain a  multiple map, pick the first (match does this automatically)
index1 = match(featureNames(eset), map[,1])

# now use this to map to the genes GRanges
index2 = match(map[index1,2], as.character(mcols(genes)$GENEID))

# remove NAs
index3 = which(!is.na(index2))
index2 = index2[index3]

# subset objects to map
genes = genes[index2,]  
neweset =  eset[index3,]

# create a summarized experiment for eset
se = SummarizedExperiment(assays=exprs(neweset),
                          rowData=genes,
                          colData=DataFrame(pData(neweset)))
se

dim(assay(se))
length(granges(se))

tss = resize(granges(se), 1)

head(tss)

sum(start(tss) <= 50*10^6 & seqnames(tss) == "chr1")

# plot differential expression across the genome

# re-order se
se = se[order(granges(se)),]
ind = se$group==1
de = rowMeans(assay(se)[,ind]) - rowMeans(assay(se)[,!ind])
chrs = unique(seqnames(se))

# chr1-4
par(mfrow=c(3,2))
for(i in c(1:4)){
  ind = which(seqnames(se) == chrs[i])
  plot(start(se)[ind], de[ind], ylim=c(-1,1), main=as.character(chrs[i]))
  abline(h=0)
}

# chrX and Y
for(i in 23:24){
  ind = which(seqnames(se) == chrs[i])
  plot(start(se)[ind], de[ind], ylim=c(-5,5), main=as.character(chrs[i]))
  abline(h=0)
}
