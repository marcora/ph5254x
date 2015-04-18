# source("http://bioconductor.org/biocLite.R")
# biocLite("AnnotationHub")

library(AnnotationHub)

ah = AnnotationHub()

ah

query(ah, "HepG2")

###############################################################################

library(org.Hs.eg.db)
library(Homo.sapiens)

org.Hs.eg.db # is a mapping of Entrez IDs to other resource IDs (e.g., ENSEMBL IDs, PUBMED IDs, etc.)
Homo.sapiens

keytypes(org.Hs.eg.db)
columns(org.Hs.eg.db)

keytypes(Homo.sapiens)
columns(Homo.sapiens)

select(org.Hs.eg.db, keytype="SYMBOL", key="ORMDL3", columns="PMID")
select(Homo.sapiens, keytype="SYMBOL", key="ORMDL3", columns="PMID")

select(org.Hs.eg.db, keytype="SYMBOL", key="ORMDL3", columns="GO")
select(Homo.sapiens, keytype="SYMBOL", key="ORMDL3", columns=c("GO", "TERM"))


# The Gene Ontology and the GO Annotation db group genes by label
library(GO.db)
GO.db

# access underlying SQLite server
con = GO_dbconn()
dbListTables(con)


# The KEGG db group genes by pathway
library(KEGGREST)
library(png)
library(grid)
brpng = keggGet("hsa05212", "image")
grid.raster(brpng)
