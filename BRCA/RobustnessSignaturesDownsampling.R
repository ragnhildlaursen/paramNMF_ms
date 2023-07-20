rm(list=ls())
setwd("~/projects/paramNMF_ms")
################################################################
## Robustness of Signatures
################################################################

source("BRCA/loadBRCA214models.R")
load("C:/Users/au543194/Documents/projects/paramNMF_ms/BRCA/result/notUseResults/BRCA214sig8/BRCA214optimal8sig.RData")

MonoRes = resFactors[[1]]
DiRes = resFactors[[2]]
TriRes = resFactors[[3]]
MixRes = resFactors[[4]]


mMono = cosMatch(TriRes$signatures, MonoRes$signatures)$match
mMix = cosMatch(TriRes$signatures, MixRes$signatures)$match
mDi = cosMatch(TriRes$signatures, DiRes$signatures)$match

##-----------------------------------------------------------------------
## Reestimating signatures and exposures
##-----------------------------------------------------------------------
sampleV = V
nG = dim(V)[1]     # Number of patients (genomes)
nMtTps = dim(V)[2] # Number of mutation types
noSig = nrow(TriRes$signatures)

init = 100   # number of initialisations
smallIter = 500
tol = 1e-8
nSim = 50

ResCosineMatTri = matrix(0,nrow=noSig,ncol=nSim) 
ResCosineMatMix = matrix(0,nrow=noSig,ncol=nSim) 
ResCosineMatDi = matrix(0,nrow=noSig,ncol=nSim) 
ResCosineMatMono = matrix(0,nrow=noSig,ncol=nSim)

for(nsim in 1:nSim){
  cat(nsim,"out of",nSim,"\n")
  
  estimate = TriRes$exposures%*%TriRes$signatures
  sampleV = matrix(rpois(nG*nMtTps, estimate), nrow = nG)
  ResultFixTri <- nmfprm(data = sampleV, noSignatures = noSig, designMatrices = rep(list(Mfull),noSig), tol=tol,initial = init, smallIter = smallIter)
  
  estimate = MixRes$exposures%*%MixRes$signatures
  sampleV = matrix(rpois(nG*nMtTps, estimate), nrow = nG)
  ResultFixMix <- nmfprm(data = sampleV, noSignatures = noSig, designMatrices = list(Mmono,Mdi,Mdi,Mdi,Mdi,Mdi,Mfull,Mfull),tol=tol,initial = init, smallIter = smallIter)
  
  estimate = DiRes$exposures%*%DiRes$signatures
  sampleV = matrix(rpois(nG*nMtTps, estimate), nrow = nG)
  ResultFixDi <- nmfprm(data = sampleV, noSignatures = noSig, designMatrices = rep(list(Mdi),noSig),tol=tol,initial = init, smallIter = smallIter)
  
  estimate = MonoRes$exposures%*%MonoRes$signatures
  sampleV = matrix(rpois(nG*nMtTps, estimate), nrow = nG)
  ResultFixMono <- nmfprm(data = sampleV, noSignatures = noSig, designMatrices = rep(list(Mmono),noSig),tol=tol,initial = init, smallIter = smallIter)
  
  # Compare signatures: Cosine similarity for each patient
  ResCosineMatTri[,nsim] = cosMatch(TriRes$signatures,ResultFixTri$signatures)$cossim
  ResCosineMatMix[,nsim] = cosMatch(MixRes$signatures[mMix,],ResultFixMix$signatures)$cossim
  ResCosineMatDi[,nsim] = cosMatch(DiRes$signatures[mDi,],ResultFixDi$signatures)$cossim
  ResCosineMatMono[,nsim] = cosMatch(MonoRes$signatures[mMono,],ResultFixMono$signatures)$cossim
  
  save(ResCosineMatTri,ResCosineMatMix, ResCosineMatDi,ResCosineMatMono, file = "SimilaritySignaturesBRCA214optimal8.RData")
}
