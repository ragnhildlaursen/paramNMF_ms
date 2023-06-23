rm(list=ls())
################################################################
## Robustness of exposures
################################################################
setwd("~/projects/paramNMF_ms")

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
  ResultFixTri <- NMFglmSQR(Data = sampleV, NoSignatures = noSig, Signatures = TriRes$signatures,tol=tol)
  ResultFixMix <- NMFglmSQR(Data = sampleV, NoSignatures = noSig, Signatures = MixRes$signatures,tol=tol)
  ResultFixDi <- NMFglmSQR(Data = sampleV, NoSignatures = noSig, Signatures = DiRes$signatures,tol=tol)
  ResultFixMono <- NMFglmSQR(Data = sampleV, NoSignatures = noSig, Signatures = MonoRes$signatures,tol=tol)
  
  ## Compare exposures: Cosine similarity for each patient
  ResCosineMatTri[,nsim] = cosMatch(t(TriRes$exposures),t(ResultFixTri$Exposures))$cossim
  ResCosineMatDi[,nsim] = cosMatch(t(DiRes$exposures[,mDi]),t(ResultFixDi$Exposures))$cossim
  ResCosineMatMix[,nsim] = cosMatch(t(MixRes$exposures[,mMix]),t(ResultFixMix$Exposures))$cossim
  ResCosineMatMono[,nsim] = cosMatch(t(MonoRes$exposures[,mMono]),t(ResultFixMono$Exposures))$cossim
}
save(ResCosineMatTri,ResCosineMatDi,ResCosineMatMono, file = paste0("ExposureBRCA214dwn",dwn[d]*100,".RData"))
}
t(TriRes$exposures)
