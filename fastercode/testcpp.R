rm(list=ls())
#########################################################
## Testing new NMF functions
#########################################################


library(NMF)
library(Rcpp)
library(RcppArmadillo)
library(bench)

setwd("~/projects/paramNMF_ms/")
source("BRCA/loadBRCAmodels.R")
sourceCpp("fastercode/NMF2.cpp")
sourceCpp("fastercode/NMF.cpp")

# Test 
start = Sys.time()
res = nmfprm(V, list(Mmono,Mdi,Mdi,Mmono), 4, initial = 10, maxiter = 5000)
#res2 = NMFglmSQR(V,4,list(Mmono,Mdi,Mdi,Mmono), initial = 10, maxIter = 5000)
end = Sys.time()
end - start

res1 = nmf2(V,4, tolerance = 1e-10 ,maxiter = 500); cat("cpp gkl: ", res1$gkl)

res2 = NMFglmSQR(V,4, maxIter = 5000, initial = 20, tolerance = 1e-10)
cat("sqr gkl: ", res2$gkl)

