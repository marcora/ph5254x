# source("http://bioconductor.org/biocLite.R")
# or
# library(BiocInstaller)
# biocLite("COPDSexualDimorphism.data")

library(tidyr)
library(dplyr)
library(magrittr)

library(COPDSexualDimorphism.data)
data(lgrc.expr.meta)

expr.meta %>%
  filter(gender == "2-Female") %>%
  count()

expr.meta %>%
  summarize(median(pkyrs)) %>%
  print

qqnorm(expr.meta$pkyrs)
qqline(expr.meta$pkyrs)

# source("http://bioconductor.org/biocLite.R")
# or
# library(BiocInstaller)
# biocLite("BSgenome.Hsapiens.UCSC.hg19")

library(BSgenome.Hsapiens.UCSC.hg19)
BSgenome.Hsapiens.UCSC.hg19

chr11seq = BSgenome.Hsapiens.UCSC.hg19[["chr11"]]
subseq(chr11seq, start=10^6, width=25)

sort(sapply(c("ATG","TGA","TAA","TAG"), function(pattern) {
  countPattern(pattern, chr11seq)
}))


chr7seq = BSgenome.Hsapiens.UCSC.hg19[["chr7"]]

counts = alphabetFrequency(chr7seq)

freqs = counts/sum(counts) # same as freqs = alphabetFrequency(chr7seq, as.prob = TRUE)

# source("http://bioconductor.org/biocLite.R")
# or
# library(BiocInstaller)
# biocLite("SNPlocs.Hsapiens.dbSNP.20120608")

library(SNPlocs.Hsapiens.dbSNP.20120608)

chr17snp = getSNPlocs("ch17")

head(chr17snp)

chr17snp %>%
  filter(RefSNP_id == "73971683") %$% loc
# same as: subset(chr17snp, RefSNP_id == "73971683")$loc
# same as: rsid2loc("rs73971683", caching=TRUE)

###############################################################################

# library(devtools)
# install_github("genomicsclass/tissuesGeneExpression")


library(tissuesGeneExpression)
data(tissuesGeneExpression)

head(e[,1:5])
table(tissue)


boxplot(e["209169_at",] ~ tissue, las=2)

ids = c("201884_at", "209169_at", "206269_at", "207437_at", "219832_s_at", "212827_at")

par(mfrow =c(2,3))

for(id in ids) {
  boxplot(e[id,] ~ tissue, las=2, main=id)
}

par(mfrow =c(1,1))

###############################################################################

library(Biobase)
# library(devtools)
# install_github("genomicsclass/GSE5859")
library(GSE5859)

data(GSE5859)

class(e) # e is an expressionSet
dim(e) # 8793 features for 208 samples
e


dat = exprs(e)
dim(dat) # 8793 rows (features) x 208 columns (samples)
head(dat)


sampleInfo = pData(e)
dim(sampleInfo)
head(sampleInfo)


# source("http://bioconductor.org/biocLite.R")
# or
# library(BiocInstaller)
# biocLite("hgfocus.db")

library(hgfocus.db)
annot = select(hgfocus.db, 
               keys=featureNames(e), 
               keytype="PROBEID", 
               columns=c("CHR", "CHRLOC", "SYMBOL"))
## here we pick one column from the annotation
annot = annot[match(featureNames(e), annot$PROBEID), ]
dim(annot)
head(annot)

