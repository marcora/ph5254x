library(tidyr)
library(dplyr)
library(magrittr)

library(Biobase)

?ExpressionSet

data(sample.ExpressionSet)

ses = sample.ExpressionSet

ses

dim(ses)

attributes(ses)

exprs(ses)[1:6, 1:6] # expression values

pData(ses)

experimentData(ses)


