<<<<<<< HEAD
library(devtools)
install_github("genomicsclass/GSE5859Subset")

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

=======
# library(devtools)
# install_github("genomicsclass/ERBS")
library(ERBS)

data(HepG2)
data(GM12878)

HepG2 # HepG2 cell line
GM12878 # GM12878 cell line

# overlaps
fo = findOverlaps(HepG2, GM12878)
erbs = HepG2[queryHits(fo)]
erbs = granges(erbs)

# intersection
erbs2 = intersect(HepG2, GM12878)
erbs2 = granges(erbs2)

# sort
erbs = sort(erbs)
erbs2 = sort(erbs2)

# confirm same chr and strand
all(seqnames(erbs) == seqnames(erbs2))
all(strand(erbs) == strand(erbs2))

# proportion of identical ranges based on start and end
mean(start(erbs) == start(erbs2) & end(erbs2) == end(erbs3))

# the intersection should be smaller
all(width(erbs2) <= width(erbs))

###############################################################################
library(Homo.sapiens)

ghs = genes(Homo.sapiens)

tss = resize(ghs, 1)

start(tss["100113402"])

fo = findOverlaps(HepG2, GM12878)
erbs = HepG2[queryHits(fo)]

index = nearest(erbs[4], tss)

geneid = names(tss[index])
geneid

keytypes(Homo.sapiens)
columns(Homo.sapiens)

select(Homo.sapiens, keytype="GENEID", keys=geneid, columns = c("SYMBOL"))
>>>>>>> f24ad0feac8ed22e87b68bbd07eec4c2e257ce55
