library(AnnotationHub)

ah = AnnotationHub()

mah = metadata(ah)

names(mah)

sort(table(mah$Species), decreasing=TRUE)[1:10]

names(query(query(ah, "HepG2"), "CTCF"))
length(names(query(query(ah, "HepG2"), "CTCF")))
