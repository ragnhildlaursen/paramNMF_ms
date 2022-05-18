rm(list=ls())
setwd("~/projects/paramNMF_ms")
################################################################
## Robustness of Signatures
################################################################

source("UCUT/loadUCUTmodels.R")

MonoRes = resFactors[[1]]
DiRes = resFactors[[19]]
PentaRes = resFactors[[21]]

mMono = cosMatch(PentaRes$signatures, MonoRes$signatures)$match
mDi = cosMatch(PentaRes$signatures, DiRes$signatures)$match

##-----------------------------------------------------------------------
## Reestimating signatures and exposures
##-----------------------------------------------------------------------
sampleV <- V
nG <- dim(V)[1]     # Number of patients (genomes)
nMtTps <- dim(V)[2] # Number of mutation types
nMt <- rowSums(V)
noSig = 2

nSim <- 50  # 5 if small study. 100 if intermediate. 200 if large.
init = 100   # number of initialisations
tol = 1e-8
sI = 100
mI = 10000

ResCosineMatPenta <- matrix(0,nrow=noSig,ncol=nSim) 
ResCosineMatDi <- matrix(0,nrow=noSig,ncol=nSim) 
ResCosineMatMono <- matrix(0,nrow=noSig,ncol=nSim)

for(nsim in 1:nSim){
  cat(nsim,"out of",nSim,"\n")
  nSimMt <- rowSums(V) # could potentially downsample
  
  estimate = PentaRes$exposures%*%PentaRes$signatures
  sampleV = matrix(rpois(nG*nMtTps, estimate), nrow = nG)
  sampleV[sampleV == 0] = .Machine$double.eps # Zero entries replaced with small epsilon to avoid division by zero in EM-algorithm
  ResultFixPenta <- nmfprm(data = sampleV, noSignatures = noSig, designMatrices = rep(list(Mfull),noSig), 
                           tol=tol,initial = init, smallIter = sI, maxiter = mI)
  
  estimate = DiRes$exposures%*%DiRes$signatures
  sampleV = matrix(rpois(nG*nMtTps, estimate), nrow = nG)
  sampleV[sampleV == 0] = .Machine$double.eps
  ResultFixDi <- nmfprm(data = sampleV, noSignatures = noSig, designMatrices = rep(list(Mnghbr),noSig),
                        tol=tol,initial = init, smallIter = sI, maxiter = mI)
  
  estimate = MonoRes$exposures%*%MonoRes$signatures
  sampleV = matrix(rpois(nG*nMtTps, estimate), nrow = nG)
  sampleV[sampleV == 0] = .Machine$double.eps
  ResultFixMono <- nmfprm(data = sampleV, noSignatures = noSig, designMatrices = rep(list(Mmono),noSig),
                          tol=tol,initial = init, smallIter = sI, maxiter = mI)
  
  # Compare signatures: Cosine similarity for each patient
  ResCosineMatPenta[,nsim] = cosMatch(PentaRes$signatures,ResultFixPenta$signatures)$cossim
  ResCosineMatDi[,nsim] = cosMatch(DiRes$signatures[mDi,],ResultFixDi$signatures)$cossim
  ResCosineMatMono[,nsim] = cosMatch(MonoRes$signatures[mMono,],ResultFixMono$signatures)$cossim
  
  save(ResCosineMatPenta, ResCosineMatDi, ResCosineMatMono, file = "SimilaritySignaturesUCUTparamboot.RData")
}


