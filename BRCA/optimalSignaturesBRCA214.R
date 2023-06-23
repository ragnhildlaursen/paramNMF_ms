rm(list=ls())
####################################
## Find optimal models 
####################################
library(Rcpp)
library(RcppArmadillo)
library(bench)
library(gtools)
setwd("~/projects/paramNMF_ms/")
source("BRCA/loadBRCA214models.R")

## List of 15 models for the 3 parametrizations with X signatures
noSignatures = 8

models = list(Mmono,Mdi,Mfull)
ModelCombinations = matrix(c(rep(1,8),rep(2,8),rep(3,8),c(1,rep(2,5),3,3)), nrow = 4, byrow = T)
nModels <- nrow(ModelCombinations)
MList = lapply(1:nModels,function(x) lapply(ModelCombinations[x,], function(y) models[[y]]))

##-----------------------------------------------------------
## Running the different models
##-----------------------------------------------------------

# matrix for results 
resMat <- matrix(0,nrow=nModels,ncol=noSignatures + 2)
colnames(resMat) <- c(paste0("nprm", c(1:noSignatures)),"nprmtot","GKL")
resFactors = list()
nprm = numeric(noSignatures)
init = 100
smallIter = 500
low.tolerance <- 1e-5


for (m in 1:nModels){
  for(sig in 1:noSignatures){
    nprm[sig] = ncol( MList[[m]][[sig]] )
  }
  resMat[m,c(1:noSignatures)] <- nprm
  resMat[m,"nprmtot"] <- sum(nprm)
  
  ## List of models chosen
  cat("Model:",m,"","\n")
  res <- nmfprm(data=V,noSignatures=noSignatures,designMatrices=MList[[m]],
                tolerance=low.tolerance, initial = init, maxiter = 10000, smallIter = smallIter)
  resFactors[[m]] = res
  
  ## print results
  resMat[m,"GKL"] <- res$gkl
  cat("Final result:","nprm:",nprm,", GKL:",res$gkl,"\n")
}
save(resFactors,resMat, file = "BRCA214optimal8sig.RData")

#######################################################################

