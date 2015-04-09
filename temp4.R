library(Biobase)
data(sample.ExpressionSet)
sample.ExpressionSet
samp = sample.ExpressionSet

dim(samp)

pData(samp)


pData(samp)$sex

samp$sex

dat = exprs(samp)

fems = samp[, samp$sex == 'Female']

dim(samp)
dim(fems)

sum(exprs(fems)[1,])

experimentData(samp)

annotation(samp)

cor(samp$score, exprs(samp)["31489_at",])
