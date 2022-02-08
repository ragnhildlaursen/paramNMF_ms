##-----------------------------------------------------
## UCUT 26 data analysis
##-----------------------------------------------------
source("GitModelSelection.R")

## Load UCUT data
load("UCUT_5_all.RData")

V5 = t(Vall)   # data (No. patients) x (No. Mutation types)
## UCUT has 1536 mutation types and 26 patients 
V5[V5 == 0] = .Machine$double.eps # Zero entries replaced with small epsilon to avoid division by zero in EM-algorithm
##--------------------------------------------------------
## Factors
##--------------------------------------------------------
## First and second left flanking nucleotide 
l1 = factor(substr(colnames(V5), start = 1, stop = 1))
l2 = factor(substr(colnames(V5), start = 9, stop = 9))
## Actual point mutation
m = factor(substr(colnames(V5), start = 3,stop = 5))
## First and second right flanking nucleotide 
r1 = factor(substr(colnames(V5), start = 7, stop = 7))
r2 = factor(substr(colnames(V5), start = 11, stop = 11))

##--------------------------------------------------------
## Parametrizations of a signature
##--------------------------------------------------------
Mmono = model.matrix(~0+l2+l1+m+r1+r2)          # Mono-nucleotide model
Mdi = model.matrix(~0+l2*m+l1*m+m*r1+m*r2)      # Di-nucleotide interaction with mutation
Mblend = model.matrix(~0+l2+l1*m+m*r1+r2)       # Blended mono-di-nucleotide model
Mcombi = model.matrix(~0+l2+l1*m*r1+r2)         # Combined mono-tri-nucleotide model
Mtri =  model.matrix(~0+l1*m*r1)                # Tri-nucleotide model (only one flanking nucleotide)
Mnghbr = model.matrix(~0+l2*l1+l1*m+m*r1+r1*r2) # Di-nucleotide interaction with neighbour
Mfull = model.matrix(~0+l2*l1*m*r1*r2)          # Full parametrized model
## List of the 21 models
MList <- list(list(Mmono,Mmono),
              list(Mdi,Mdi),list(Mmono,Mdi),
              list(Mblend,Mblend),list(Mmono,Mblend),list(Mdi,Mblend),
              list(Mcombi,Mcombi),list(Mmono,Mcombi),list(Mdi,Mcombi),list(Mblend,Mcombi),
              list(Mtri,Mtri),list(Mmono,Mtri),list(Mdi,Mtri),list(Mblend,Mtri),
              list(Mcombi,Mtri),
              list(Mnghbr,Mnghbr),list(Mmono,Mnghbr),list(Mdi,Mnghbr),list(Mblend,Mnghbr),
              list(Mcombi,Mnghbr),list(Mtri,Mnghbr))
nModels <- length(MList)
##-----------------------------------------------------------
## Matrix with summary of the results
##-----------------------------------------------------------
resMat <- matrix(0,nrow=nModels,ncol=5)
colnames(resMat) <- c("nprm1","nprm2","nprmtot","seed","GKL")
## First run EM algorithm 30 times with high tolerance to 
## identify global minimum 
## Second run EM algorithm with seed from global minimum but 
## with low tolerance (i.e. until convergence)
## Here we use high tolerance=25 and low tolerance=5
## to run the algorithm in a few minutes.
## In order to obtain the results in the paper you *must* run with
## high tolerance=5 and low tolerance=1.
## On a normal laptop this takes around 12 hours. 
## The GKL is summarized in resMat. 
## We provide resMat for our run in file "GitUCUT26results.txt"
high.tolerance <- 25
low.tolerance <- 5
for (i in 1:nModels){
  nprm1 <- ncol( MList[[i]][[1]] )
  nprm2 <- ncol( MList[[i]][[2]] )
  resMat[i,1] <- nprm1
  resMat[i,2] <- nprm2
  resMat[i,3] <- nprm1+nprm2
  tmp <- rep(0,30)
  cat("Model:",i,"","\n")
  cat("EM Run: ")
  for (j in 1:length(tmp)){
    cat(j,"")
    tmp[j] <- 
      NMFglmSQR(Data=V5,DesignMatrix=MList[[i]],tolerance=high.tolerance,Seeds=j)$gkl
    ## Run the EM algorithm many times with a high tolerance to identify global minimum
  }
  ## Choose the best initial value and run again until convergence, i.e.
  ## stop when the tolerance is small
  cat("\n")
  jmin <- which.min(tmp)
  resMat[i,4] <- jmin
  res <- NMFglmSQR(Data=V5,DesignMatrix=MList[[i]],tolerance=low.tolerance,Seeds=jmin)$gkl
  ## 
  resMat[i,5] <- res
  cat("Final result:","nprm1:",nprm1,"; nprm2:",nprm2,
      "; nprmtotal:",nprm1+nprm2,"; seed:",jmin,", GKL:",res,"\n" )
}
resMat <- resMat[sort(resMat[,"nprmtot"],index=TRUE)$ix,]
print(resMat)

