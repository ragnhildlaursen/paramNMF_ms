rm(list=ls())
################################################################
## Robustness of exposures
################################################################
setwd("~/projects/paramNMF_ms")
source("UCUT/loadUCUTmodels.R")

MonoRes = resFactors[[1]]
DiRes = resFactors[[19]]
PentaRes = resFactors[[21]]

mMono = cosMatch(PentaRes$signatures, MonoRes$signatures)$match
mDi = cosMatch(PentaRes$signatures, DiRes$signatures)$match

##-----------------------------------------------------------------------
## Sample a number of mutations from each patient (downsampling)
##-----------------------------------------------------------------------
sampleV <- V
nG <- dim(V)[1]     # Number of patients (genomes)
noSig = 2
nMtTps <- dim(V)[2] # Number of mutation types
nMt <- rowSums(V)
dwn <- c(0.01,0.02,0.05)  # Down sampling 100, 20, 10
nSim <- 50  # 5 if small study. 100 if intermediate. 200 if large.
tol = 1e-2

for(d in 1:length(dwn)){
ResCosineMatPenta <- matrix(0,nrow=noSig,ncol=nSim) 
ResCosineMatDi <- matrix(0,nrow=noSig,ncol=nSim) 
ResCosineMatMono <- matrix(0,nrow=noSig,ncol=nSim)

for(nsim in 1:nSim){
  cat(nsim,"out of",nSim,"\n")
  nSimMt <- round(rowSums(V)*dwn[d])
  for (i in 1:nG){
    set.seed(NULL)
    sim <- sample( x=1:nMtTps , size=nSimMt[i], replace=TRUE, prob=V[i,])
    sampleV[i,] <- tabulate(sim,nbins=nMtTps)
  }
  
  ResultFixPenta <- NMFglmSQR(Data = sampleV, NoSignatures = noSig, Signatures = PentaRes$signatures,tol=tol)
  ResultFixDi <- NMFglmSQR(Data = sampleV, NoSignatures = noSig, Signatures = DiRes$signatures,tol=tol)
  ResultFixMono <- NMFglmSQR(Data = sampleV, NoSignatures = noSig, Signatures = MonoRes$signatures,tol=tol)
  
  ## Compare exposures: Cosine similarity for each patient
  ResCosineMatPenta[,nsim] = cosMatch(t(PentaRes$exposures),t(ResultFixPenta$Exposures))$cossim
  ResCosineMatDi[,nsim] = cosMatch(t(DiRes$exposures[,mDi]),t(ResultFixDi$Exposures))$cossim
  ResCosineMatMono[,nsim] = cosMatch(t(MonoRes$exposures[,mMono]),t(ResultFixMono$Exposures))$cossim
}
save(ResCosineMatPenta,ResCosineMatDi,ResCosineMatMono, file = paste0("ExposureUCUTdwn",dwn[d]*100,".RData"))
}


