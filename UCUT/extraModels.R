############################################
## Extra models checked
############################################

setwd("~/projects/paramNMF_ms/")
source("UCUT/loadUCUTmodels.R")







init = 100
smallIter = 500
low.tolerance <- 1e-8

##-----------------------------------------------------------
## Di- and Tri- nucleotide model (L2*L1 + L1*M*R1 + R1*R2)
##-----------------------------------------------------------

# start = Sys.time()
# resditri = nmfprm(data=V,noSignatures=2,designMatrices=list(Mditri,Mditri),
#                tolerance=low.tolerance, initial = init, maxiter = 10000, smallIter = smallIter)
# end = Sys.time()
# end - start

## Together with other models ----------------------------------------
noSignatures = 2
nModels = 6
MList = list(list(Mditri,Mmono),
             list(Mditri,Mblend),
             list(Mditri,Mcombi),
             list(Mditri,Mtri),
             list(Mditri,Mnghbr),
             list(Mditri,Mfull))


resMat <- matrix(0,nrow=nModels,ncol=4)
colnames(resMat) <- c("nprm1","nprm2","nprmtot","GKL")
resFactors = list()

init = 500
smallIter = 500
low.tolerance <- 1e-8

for(m in 1:nModels){
  nprm1 <- ncol( MList[[m]][[1]] )
  nprm2 <- ncol( MList[[m]][[2]] )
  resMat[m,1] <- nprm1
  resMat[m,2] <- nprm2
  resMat[m,3] <- nprm1+nprm2
  
  cat("Model:",m,"out of",nModels,"\n")
  res = nmfprm(data=V,noSignatures=2,designMatrices=MList[[m]],
               tolerance=low.tolerance, initial = init, maxiter = 10000, smallIter = smallIter)
  
  resFactors[[m]] = res
  
  ## 
  resMat[m,4] <- res$gkl
  cat("Final result:","nprm1:",nprm1,"; nprm2:",nprm2,
      "; nprmtotal:",nprm1+nprm2,", GKL:",res$gkl,"\n")
  
  save(resMat,resFactors, file = "UCUTditriMIXmodel.RData")
}



# load("result/UCUTditrimodel.RData")
# dim(Mditri)
# 2*resditri$gkl + log(5260)*240
# 2*9008 + 1131