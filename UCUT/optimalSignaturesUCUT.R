#rm(list=ls())
################################################################
## Find optimal signatures for UCUT
################################################################

setwd("~/projects/paramNMF_ms/")
source("UCUT/loadUCUTmodels.R")

##-----------------------------------------------------------
## Matrix with summary of the results
##-----------------------------------------------------------
resMat <- matrix(0,nrow=nModels,ncol=4)
colnames(resMat) <- c("nprm1","nprm2","nprmtot","GKL")
resFactors = list()

low.tolerance <- 1e-8

for(m in 1:nModels){
  nprm1 <- ncol( MList[[m]][[1]] )
  nprm2 <- ncol( MList[[m]][[2]] )
  resMat[m,1] <- nprm1
  resMat[m,2] <- nprm2
  resMat[m,3] <- nprm1+nprm2
  
  cat("Model:",m,"out of",nModels,"\n")
  res = nmfprm(data=V,noSignatures=2,designMatrices=MList[[m]],
         tolerance=low.tolerance, initial = 500, maxiter = 10000, smallIter = 500)
  
  resFactors[[m]] = res

  ## 
  resMat[m,4] <- res$gkl
  cat("Final result:","nprm1:",nprm1,"; nprm2:",nprm2,
      "; nprmtotal:",nprm1+nprm2,", GKL:",res$gkl,"\n")
}

save(resMat, file = "UCUTmodelsummary.RData")
save(resFactors, file = "UCUTmodelFactors.RData")


