setwd("~/projects/paramNMF_ms/")
rm(list=ls())
####################################
## Find optimal models BRCA
####################################
library(Rcpp)
library(RcppArmadillo)
library(bench)
library(gtools)

library(foreach)
library(doParallel)
library(parallel)
library(SQUAREM)
library(MASS)

#cores=detectCores()
#cl <- makeCluster(cores[1]-1)

source("UCUT/loadUCUTmodels.R")

init = 1
smallIter = 100
low.tolerance <- 1e-1
maxiter = 1000
sig = c(1:3)
models = list(Mmono,Mtri,Mblend,Mcombi,Mnghbr,Mditri)

mcomb = matrix(0, nrow = length(sig), ncol = 2)
colnames(mcomb) = c("sig","nomodel")
mcomb[,1] = sig

modellist = list()
## create list of model combinations
for(i in 1:length(sig)){
  noSignatures = mcomb[i,1]
  
  modellist[[i]] = combinations(length(models), noSignatures, 
                                repeats.allowed = TRUE)
  mcomb[i,2] = nrow(modellist[[i]])
}

#registerDoParallel(cl)

allmodels = foreach(s=1:length(sig), .packages=c('Rcpp','RcppArmadillo'), .export = ls(globalenv())) %:%
  foreach(m=1:mcomb[s,2], .packages=c('Rcpp','RcppArmadillo'), .export = ls(globalenv()), .combine='cbind') %do% {
    sourceCpp("NMF2.cpp")
    MList = lapply(modellist[[s]][m,], function(y) models[[y]])
    
    noSignatures = sig[s]
    nprm = numeric(noSignatures)
    for(i in 1:noSignatures){
      nprm[i] = ncol( MList[[i]] )
    }
    
    ## List of models chosen
    res = nmfprm(data=V,noSignatures=noSignatures,designMatrices=MList,
                 tolerance=low.tolerance, initial = init, maxiter = maxiter, smallIter = smallIter)
    
    bic = res$gkl*2 + (nrow(V)*sig[s] + sum(nprm))*log(sum(V > 0))
    
    ## get results
    c(noSignatures,nprm, sum(nprm),res$gkl,bic)
  }
#stopCluster(cl)

save(allmodels, file = "UCUTallmodels1init.RData")