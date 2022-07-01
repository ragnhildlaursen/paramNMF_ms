rm(list=ls())
####################################
## Find optimal models 
####################################
library(Rcpp)
library(RcppArmadillo)
library(bench)
library(gtools)
setwd("~/projects/paramNMF_ms/")
source("BRCA/loadBRCA21models.R")

## List of 15 models for the 3 parametrizations with X signatures
noSignatures = 4

models = list(Mmono,Mdi,Mfull)
ModelCombinations = combinations(length(models), noSignatures, 
                                 repeats.allowed = TRUE)
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
init = 500
smallIter = 500
low.tolerance <- 1e-8

 
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


#######################################################################

