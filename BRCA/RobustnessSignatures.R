rm(list=ls())
setwd("~/projects/paramNMF_ms")
################################################################
## Robustness of Signatures
################################################################
# library(foreach)
# library(doParallel)
# library(parallel)
# 
# cores=detectCores()
# cl <- makeCluster(cores[1]-1)

source("BRCA/loadBRCA214models.R")

MonoRes = resFactors[[1]]
MixRes = resFactors[[13]]
DiRes = resFactors[[37]]
TriRes = resFactors[[45]]

mMono = cosMatch(TriRes$signatures, MonoRes$signatures)$match
mMix = cosMatch(TriRes$signatures, MixRes$signatures)$match
mDi = cosMatch(TriRes$signatures, DiRes$signatures)$match

##-----------------------------------------------------------------------
## Sample a number of mutations from each patient (downsampling)
##-----------------------------------------------------------------------
sampleV <- V
nG <- dim(V)[1]     # Number of patients (genomes)
nMtTps <- dim(V)[2] # Number of mutation types
nMt <- rowSums(V)
noSig = 8

nSim <- 50  # 5 if small study. 100 if intermediate. 200 if large.
init = 50   # number of initialisations
tol = 1e-6

# TriRes <- NMFglmSQR(Data = sampleV, NoSignatures = noSig, DesignMatrix = rep(list(Mfull),noSig), tolerance = tol,initial = init)
# MixRes <- NMFglmSQR(Data = sampleV, NoSignatures = noSig, DesignMatrix = list(Mfull,Mdi,Mdi,Mmono),tol=tol,initial = init)
# DiRes <- NMFglmSQR(Data = sampleV, NoSignatures = noSig, DesignMatrix = rep(list(Mdi),noSig),tolerance = tol,initial = init)
# MonoRes <- NMFglmSQR(Data = sampleV, NoSignatures = noSig, DesignMatrix = rep(list(Mmono),noSig),tolerance = tol,initial = init)

#save(TriRes,DiRes,MonoRes, file = "FactorsBRCA119.RData")
#cat("GKL Di:", DiRes$gkl, " with convergence ", DiRes$conv, "\n")
#cat("GKL Mix:", MixRes$gkl, " with convergence ", MixRes$conv, "\n")

ResCosineMatTri <- matrix(0,nrow=noSig,ncol=nSim) 
ResCosineMatMix <- matrix(0,nrow=noSig,ncol=nSim) 
ResCosineMatDi <- matrix(0,nrow=noSig,ncol=nSim) 
ResCosineMatMono <- matrix(0,nrow=noSig,ncol=nSim)

for(nsim in 1:nSim){
  cat(nsim,"out of",nSim,"\n")
  nSimMt <- rowSums(V) # could potentially downsample
  for (i in 1:nG){
    set.seed(NULL)
    sim <- sample( x=1:nMtTps , size=nSimMt[i], replace=TRUE, prob=V[i,])
    sampleV[i,] <- tabulate(sim,nbins=nMtTps)
  }
  
  # ResultFixTri <- NMFglmSQR(Data = sampleV, NoSignatures = noSig, DesignMatrix = rep(list(Mfull),noSig), tol=tol,initial = init)
  # ResultFixMix <- NMFglmSQR(Data = sampleV, NoSignatures = noSig, DesignMatrix = list(Mfull,Mdi,Mdi,Mmono),tol=tol,initial = init)
  # ResultFixDi <- NMFglmSQR(Data = sampleV, NoSignatures = noSig, DesignMatrix = rep(list(Mdi),noSig),tol=tol,initial = init)
  # ResultFixMono <- NMFglmSQR(Data = sampleV, NoSignatures = noSig, DesignMatrix = rep(list(Mmono),noSig),tol=tol,initial = init)
  # cat("GKL Di:", ResultFixDi$gkl, " with convergence ", ResultFixDi$conv, "\n")
  # cat("GKL Mix:", ResultFixMix$gkl, " with convergence ", ResultFixMix$conv, "\n")
  
  ResultFixTri <- nmfprm(data = sampleV, noSignatures = noSig, designMatrices = rep(list(Mfull),noSig), tol=tol,initial = init)
  ResultFixMix <- nmfprm(data = sampleV, noSignatures = noSig, designMatrices = list(Mmono,Mmono,Mmono,Mmono,Mfull,Mfull,Mdi,Mdi),tol=tol,initial = init)
  ResultFixDi <- nmfprm(data = sampleV, noSignatures = noSig, designMatrices = rep(list(Mdi),noSig),tol=tol,initial = init)
  ResultFixMono <- nmfprm(data = sampleV, noSignatures = noSig, designMatrices = rep(list(Mmono),noSig),tol=tol,initial = init)
  # Compare signatures: Cosine similarity for each patient
  ResCosineMatTri[,nsim] = cosMatch(TriRes$signatures,ResultFixTri$signatures)$cossim
  ResCosineMatMix[,nsim] = cosMatch(MixRes$signatures[mMix,],ResultFixMix$signatures)$cossim
  ResCosineMatDi[,nsim] = cosMatch(DiRes$signatures[mDi,],ResultFixDi$signatures)$cossim
  ResCosineMatMono[,nsim] = cosMatch(MonoRes$signatures[mMono,],ResultFixMono$signatures)$cossim
  
  save(ResCosineMatTri,ResCosineMatMix, ResCosineMatDi,ResCosineMatMono, file = "SimilaritySignaturesBRCA214.RData")
}

