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
setwd("~/projects/paramNMF_ms")
library(ggplot2)
library(ggpubr)
source("BRCA/loadBRCAmodels.R")

##-----------------------------------------------------------------------
## Sample a number of mutations from each patient (downsampling)
##-----------------------------------------------------------------------
sampleV <- V
nG <- dim(V)[1]     # Number of patients (genomes)
nMtTps <- dim(V)[2] # Number of mutation types
nMt <- rowSums(V)
noSig = 4

nSim <- 10  # 5 if small study. 100 if intermediate. 200 if large.
init = 20   # number of initialisations
tolhigh = 0.2
tol = 0.1

TriRes <- NMFglmSQR(Data = V, NoSignatures = noSig, DesignMatrix = rep(list(Mfull),noSig), tol=tol,Seeds=c(1:init))
DiRes <- NMFglmSQR(Data = V, NoSignatures = noSig, DesignMatrix = rep(list(Mdi),noSig),tol=tol,Seeds=sample(1:1000,init))
MixRes <- NMFglmSQR(Data = V, NoSignatures = noSig, DesignMatrix = list(Mmono,Mdi,Mdi,Mfull),tol=tol,Seeds=sample(1:1000,init))
MonoRes <- NMFglmSQR(Data = V, NoSignatures = noSig, DesignMatrix = rep(list(Mmono),noSig),tol=tol,Seeds=c(1:init))

ResCosineMatTri <- matrix(0,nrow=noSig,ncol=nSim) 
ResCosineMatDi <- matrix(0,nrow=noSig,ncol=nSim) 
ResCosineMatMono <- matrix(0,nrow=noSig,ncol=nSim)
#foreach(nsim=1:nSim, .combine=c, .inorder=FALSE,.export = ls(globalenv())) %dopar% {
for(nsim in 1:nSim){
  cat(nsim,"out of",nSim,"\n")
  nSimMt <- rowSums(V)
  set.seed(nsim + 123456)
  # for (i in 1:nG){
  #   sim <- sample(x=1:nMtTps , size=nSimMt[i], replace=TRUE, prob=V[i,])
  #   sampleV[i,] <- tabulate(sim,nbins=nMtTps)
  # }
  seed = sample(1:1000,init)
  tmp = matrix(0, nrow = init, ncol = 3)
  for(s in 1:init){
  tmp[s,1] <- NMFglmSQR(Data = sampleV, NoSignatures = noSig, DesignMatrix = rep(list(Mfull),noSig), tol=tolhigh,Seeds=seed[s])$gkl
  tmp[s,2] <- NMFglmSQR(Data = sampleV, NoSignatures = noSig, DesignMatrix = rep(list(Mdi),noSig),tol=tolhigh,Seeds=seed[s])$gkl
  tmp[s,3] <- NMFglmSQR(Data = sampleV, NoSignatures = noSig, DesignMatrix = rep(list(Mmono),noSig),tol=tolhigh,Seeds=seed[s])$gkl
  }
  
  ResultFixTri <- NMFglmSQR(Data = sampleV, NoSignatures = noSig, DesignMatrix = rep(list(Mfull),noSig), tol=tol,Seeds=seed[which.min(tmp[,1])])
   
  ResultFixDi <- NMFglmSQR(Data = sampleV, NoSignatures = noSig, DesignMatrix = rep(list(Mdi),noSig),tol=tol,Seeds=seed[which.min(tmp[,2])])
  
  ResultFixMono <- NMFglmSQR(Data = sampleV, NoSignatures = noSig, DesignMatrix = rep(list(Mmono),noSig),tol=tol,Seeds=seed[which.min(tmp[,3])])
  
  ## Compare signatures: Cosine similarity for each patient
  ResCosineMatTri[,nsim] = cosMatch(TriRes$Signatures,ResultFixTri$Signatures)$cossim
  res = cosMatch(DiRes$Signatures,ResultFixDi$Signatures)
  ResCosineMatDi[,nsim] = res$cossim
  ResCosineMatMono[,nsim] = cosMatch(MonoRes$Signatures,ResultFixMono$Signatures)$cossim
  
  save(ResCosineMatTri,ResCosineMatDi,ResCosineMatMono, file = "SimilaritySignaturesBRCA2.RData")
}
#stopCluster(cl)
# load("SimilaritySignaturesBRCA.RData")
# 
# 
# Transforming data

datTri = data.frame(s = t(ResCosineMatTri), type = "tri")
m = cosMatch(TriRes$Signatures, DiRes$Signatures)$match
datDi = data.frame(s = t(ResCosineMatDi[m,]), type = "di")
m = cosMatch(TriRes$Signatures, MonoRes$Signatures)$match
datMono = data.frame(s = t(ResCosineMatMono[m,]), type = "mono")
dat1 = rbind(datTri,datDi,datMono)
dat2 = reshape(dat1, varying = colnames(dat1)[1:noSig], direction = "long")
dat2$type = factor(dat2$type, levels = c("mono","di","mix","tri"))
# remove mix
#dat2 = dat2[dat2$type != "di",]

##------------------------------------------------------------------
## Compare true exposures and estimated exposures
##------------------------------------------------------------------

p1 = ggplot(dat2[dat2$time == 1,], aes(x=type, y=s, fill=type)) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.9) +
  geom_violin(position=position_dodge(1)) +
  ggtitle("Signature 1")
p2 = ggplot(dat2[dat2$time == 2,], aes(x=type, y=s, fill=type)) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.9) +
  geom_violin(position=position_dodge(1)) +
  ggtitle("Signature 2")
p3 = ggplot(dat2[dat2$time == 3,], aes(x=type, y=s, fill=type)) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.9) +
  geom_violin(position=position_dodge(1)) +
  ggtitle("Signature 3")
p4 = ggplot(dat2[dat2$time == 4,], aes(x=type, y=s, fill=type)) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.9) +
  geom_violin(position=position_dodge(1)) +
  ggtitle("Signature 4")
ggarrange(p1,p2,p3,p4, nrow = 2, ncol = 2, common.legend = TRUE)

par(mfrow = c(4,2))
for(i in 1:4){
  plot(DiRes$Signatures[i,], type = "h")
  plot(ResultFixDi$Signatures[res$match[i],], type = "h")
}
m = cosMatch(MonoRes$Signatures,ResultFixMono$Signatures)$match
for(i in 1:4){
  plot(MonoRes$Signatures[i,], type = "h")
  plot(ResultFixMono$Signatures[i,], type = "h")
}
