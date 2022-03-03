rm(list=ls())
#setwd("C:/Users/au543194/Documents/projects/paramNMF_ms")
################################################################
## Robustness of Signatures
################################################################
# library(foreach)
# library(doParallel)
# library(parallel)
# 
# cores=detectCores()
# cl <- makeCluster(cores[1]-1)
library(ggplot2)
library(ggpubr)
source("loadBRCAmodels.R")

##-----------------------------------------------------------------------
## Sample a number of mutations from each patient (downsampling)
##-----------------------------------------------------------------------
sampleV <- V
nG <- dim(V)[1]     # Number of patients (genomes)
nMtTps <- dim(V)[2] # Number of mutation types
nMt <- rowSums(V)
noSig = 4

nSim <- 10  # 5 if small study. 100 if intermediate. 200 if large.
init = 50   # number of initialisations
tol = 0.1

TriRes <- NMFglmSQR(Data = sampleV, NoSignatures = noSig, DesignMatrix = rep(list(Mfull),noSig), tolerance = tol,initial = init)
MixRes <- NMFglmSQR(Data = sampleV, NoSignatures = noSig, DesignMatrix = list(Mfull,Mdi,Mdi,Mmono),tol=tol,initial = init)
DiRes <- NMFglmSQR(Data = sampleV, NoSignatures = noSig, DesignMatrix = rep(list(Mdi),noSig),tolerance = tol,initial = init)
MonoRes <- NMFglmSQR(Data = sampleV, NoSignatures = noSig, DesignMatrix = rep(list(Mmono),noSig),tolerance = tol,initial = init)
save(TriRes,MixRes,DiRes,MonoRes, file = "OneResultBRCA.RData")
cat("GKL Di:", DiRes$gkl, " with convergence ", DiRes$conv, "\n")
cat("GKL Mix:", MixRes$gkl, " with convergence ", MixRes$conv, "\n")

ResCosineMatTri <- matrix(0,nrow=noSig,ncol=nSim) 
ResCosineMatMix <- matrix(0,nrow=noSig,ncol=nSim) 
ResCosineMatDi <- matrix(0,nrow=noSig,ncol=nSim) 
ResCosineMatMono <- matrix(0,nrow=noSig,ncol=nSim)

#foreach(nsim=1:nSim, .combine=c, .inorder=FALSE,.export = ls(globalenv())) %dopar% {
for(nsim in 1:nSim){
  cat(nsim,"out of",nSim,"\n")
  nSimMt <- rowSums(V) # could potentially downsample
  # for (i in 1:nG){
  #   set.seed(NULL)
  #   sim <- sample( x=1:nMtTps , size=nSimMt[i], replace=TRUE, prob=V[i,])
  #   sampleV[i,] <- tabulate(sim,nbins=nMtTps)
  # }
  
  ResultFixTri <- NMFglmSQR(Data = sampleV, NoSignatures = noSig, DesignMatrix = rep(list(Mfull),noSig), tol=tol,initial = init)
  ResultFixMix <- NMFglmSQR(Data = sampleV, NoSignatures = noSig, DesignMatrix = list(Mfull,Mdi,Mdi,Mmono),tol=tol,initial = init)
  ResultFixDi <- NMFglmSQR(Data = sampleV, NoSignatures = noSig, DesignMatrix = rep(list(Mdi),noSig),tol=tol,initial = init)
  ResultFixMono <- NMFglmSQR(Data = sampleV, NoSignatures = noSig, DesignMatrix = rep(list(Mmono),noSig),tol=tol,initial = init)
  cat("GKL Di:", ResultFixDi$gkl, " with convergence ", ResultFixDi$conv, "\n")
  cat("GKL Mix:", ResultFixMix$gkl, " with convergence ", ResultFixMix$conv, "\n")
  
  # Compare signatures: Cosine similarity for each patient
  ResCosineMatTri[,nsim] = cosMatch(TriRes$Signatures,ResultFixTri$Signatures)$cossim
  ResCosineMatMix[,nsim] = cosMatch(MixRes$Signatures,ResultFixMix$Signatures)$cossim
  ResCosineMatDi[,nsim] = cosMatch(DiRes$Signatures,ResultFixDi$Signatures)$cossim
  ResCosineMatMono[,nsim] = cosMatch(MonoRes$Signatures,ResultFixMono$Signatures)$cossim
  
  save(ResCosineMatTri,ResCosineMatMix, ResCosineMatDi,ResCosineMatMono, file = "SimilaritySignaturesBRCA.RData")
}
#stopCluster(cl)
#load("SimilaritySignaturesBRCA.RData")


# # Transforming data
# 
# datTri = data.frame(s = t(ResCosineMatTri), type = "tri")
# #m = cosMatch(TriRes$Signatures, MixRes$Signatures)$match
# #datMix = data.frame(s = t(ResCosineMatMix[m,]), type = "mix")
# m = cosMatch(TriRes$Signatures, DiRes$Signatures)$match
# datDi = data.frame(s = t(ResCosineMatDi[m,]), type = "di")
# m = cosMatch(TriRes$Signatures, MonoRes$Signatures)$match
# datMono = data.frame(s = t(ResCosineMatMono[m,]), type = "mono")
# dat1 = rbind(datTri,datDi,datMono)
# dat2 = reshape(dat1, varying = colnames(dat1)[1:4], direction = "long")
# dat2$type = factor(dat2$type, levels = c("mono","di","mix","tri"))
# # remove mix 
# #dat2 = dat2[dat2$type != "mix",]
# 
# ##------------------------------------------------------------------
# ## Compare true exposures and estimated exposures
# ##------------------------------------------------------------------
# par(mfrow = c(2,2))
# p1 = ggplot(dat2[dat2$time == 1,], aes(x=type, y=s, fill=type)) +
#   geom_violin(position=position_dodge(1)) + 
#   geom_dotplot(binaxis='y', stackdir='center', dotsize=0.9) +
#   ggtitle("Signature 1")
# p2 = ggplot(dat2[dat2$time == 2,], aes(x=type, y=s, fill=type)) +
#   geom_violin(position=position_dodge(1)) + 
#   geom_dotplot(binaxis='y', stackdir='center', dotsize=0.9) +
#   ggtitle("Signature 2")
# p3 = ggplot(dat2[dat2$time == 3,], aes(x=type, y=s, fill=type)) +
#   geom_violin(position=position_dodge(1)) + 
#   geom_dotplot(binaxis='y', stackdir='center', dotsize=0.9) +
#   ggtitle("Signature 3")
# p4 = ggplot(dat2[dat2$time == 4,], aes(x=type, y=s, fill=type)) +
#   geom_violin(position=position_dodge(1)) + 
#   geom_dotplot(binaxis='y', stackdir='center', dotsize=0.9) +
#   ggtitle("Signature 4")
# ggarrange(p1,p2,p3,p4, nrow = 2, ncol = 2, common.legend = TRUE)
# 

