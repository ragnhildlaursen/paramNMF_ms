rm(list=ls())
####################################
## Find optimal models 
####################################
library(NMF)
library(Rcpp)
library(RcppArmadillo)
library(bench)
library(gtools)
setwd("~/projects/paramNMF_ms/")
source("BRCA/loadBRCA214models.R")

# find optimal number of signatures
sig = c(2:7)
idx = 1
gkl = c()
for(i in sig){
  gkl[idx] = nmfprm(V, rep(list(Mfull),i),i, initial = 20, maxiter = 10000)$gkl
  idx = idx + 1
}

BIC = 2*gkl + log(nrow(V)*ncol(V))*(nrow(V) + ncol(V))*sig
plot(sig,BIC)
which.min(BIC)
nmfprm(V, rep(list(Mdi),4),4, initial = 20, maxiter = 5000)$gkl
## List of 15 models for the 3 parametrizations with X signatures
noSignatures = 8 
range(V)

models = list(Mmono,Mdi,Mfull)
ModelCombinations = combinations(length(models), noSignatures, 
                                 repeats.allowed = TRUE)
nModels <- nrow(ModelCombinations)
MList = lapply(1:nModels,function(x) lapply(ModelCombinations[x,], function(y) models[[y]]))

##-----------------------------------------------------------
## Matrix with summary of the results
##-----------------------------------------------------------
resMat <- matrix(0,nrow=nModels,ncol=noSignatures + 2)
colnames(resMat) <- c(paste0("nprm", c(1:noSignatures)),"nprmtot","GKL")
resFactors = list()
nprm = numeric(noSignatures)

low.tolerance <- 1e-8

for (m in 1:nModels){
  for(sig in 1:noSignatures){
    nprm[sig] = ncol( MList[[m]][[sig]] )
  }
  resMat[m,c(1:noSignatures)] <- nprm
  resMat[m,"nprmtot"] <- sum(nprm)
  ## List of models chosen
  cat("Model:",m,"","\n")
  res <- nmfprm(data=V,noSignatures=8,designMatrices=MList[[m]],
                   tolerance=low.tolerance, initial = 50, maxiter = 10000)
  resFactors[[m]] = res
  
  ## print results
  resMat[m,"GKL"] <- res$gkl
  cat("Final result:","nprm:",nprm,", GKL:",res$gkl,"\n" )
}

save(resMat, file = "BRCA214modelsummary2.RData")
save(resFactors, file = "BRCA214modelFactors2.RData")
