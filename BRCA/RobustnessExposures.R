rm(list=ls())
################################################################
## Robustness of exposures
################################################################
setwd("~/projects/paramNMF_ms")
source("BRCA/loadBRCA21models.R")

MonoRes = resFactors[[1]]
MixRes = resFactors[[8]]
DiRes = resFactors[[11]]
TriRes = resFactors[[15]]

mMono = cosMatch(TriRes$signatures, MonoRes$signatures)$match
mMix = cosMatch(TriRes$signatures, MixRes$signatures)$match
mDi = cosMatch(TriRes$signatures, DiRes$signatures)$match

##-----------------------------------------------------------------------
## Reestimating exposures after downsamling the counts
##-----------------------------------------------------------------------
sampleV <- V
nG <- dim(V)[1]                 # Number of patients (genomes)
noSig = nrow(TriRes$signatures) # Number of signatures                
nMtTps <- dim(V)[2]             # Number of mutation types
nMt <- rowSums(V)               # Total number of mutations for each patients
dwn <- c(0.01,0.02,0.05)        # Down sampling percentages
nSim <- 50                      # Number of simulations
tol = 1e-2                      # tolerance for stopping



for(d in 1:length(dwn)){
ResCosineMatTri <- matrix(0,nrow=noSig,ncol=nSim) 
ResCosineMatMix <- matrix(0,nrow=noSig,ncol=nSim) 
ResCosineMatDi <- matrix(0,nrow=noSig,ncol=nSim) 
ResCosineMatMono <- matrix(0,nrow=noSig,ncol=nSim)
for (nsim in 1:nSim){
  cat(nsim,"out of",nSim,"\n")
  
  # downsampling
  nSimMt <- round(rowSums(V)*dwn[d])
  for (i in 1:nG){
    set.seed(NULL)
    sim <- sample( x=1:nMtTps , size=nSimMt[i], replace=TRUE, prob=V[i,])
    sampleV[i,] <- tabulate(sim,nbins=nMtTps)
  }
  # reestimating exposures
  ResultFixTri <- NMFglmSQR(Data = sampleV, NoSignatures = noSig, Signatures = TriRes$Signatures,tol=tol)
  ResultFixMix <- NMFglmSQR(Data = sampleV, NoSignatures = noSig, Signatures = MixRes$Signatures,tol=tol)
  ResultFixDi <- NMFglmSQR(Data = sampleV, NoSignatures = noSig, Signatures = DiRes$Signatures,tol=tol)
  ResultFixMono <- NMFglmSQR(Data = sampleV, NoSignatures = noSig, Signatures = MonoRes$Signatures,tol=tol)
  
  ## Compare exposures: Cosine similarity for each patient
  ResCosineMatTri[,nsim] = cosMatch(t(TriRes$Exposures),t(ResultFixTri$Exposures))$cossim
  ResCosineMatDi[,nsim] = cosMatch(t(DiRes$Exposures[,mDi]),t(ResultFixDi$Exposures))$cossim
  ResCosineMatMix[,nsim] = cosMatch(t(MixRes$Exposures[,mMix]),t(ResultFixMix$Exposures))$cossim
  ResCosineMatMono[,nsim] = cosMatch(t(MonoRes$Exposures[,mMono]),t(ResultFixMono$Exposures))$cossim
}
save(ResCosineMatTri,ResCosineMatDi,ResCosineMatMono, file = paste0("ExposureBRCAdwn",dwn[d]*100,".RData"))
}
