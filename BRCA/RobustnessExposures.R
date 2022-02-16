rm(list=ls())
################################################################
## Robustness of exposures
################################################################

source("loadBRCAmodels.R")
##-----------------------------------------------------------------------
## Sample a number of mutations from each patient (downsampling)
##-----------------------------------------------------------------------
sampleV <- V
nG <- dim(V)[1]     # Number of patients (genomes)
nMtTps <- dim(V)[2] # Number of mutation types
nMt <- rowSums(V)
dwn <- 20  # Down sampling 100, 50, 20, 10
nSim <- 10  # 5 if small study. 100 if intermediate. 200 if large.
ResCosineMatTri <- matrix(0,nrow=nG,ncol=nSim) 
ResCosineMatMix <- matrix(0,nrow=nG,ncol=nSim) 
ResCosineMatDi <- matrix(0,nrow=nG,ncol=nSim) 
ResCosineMatMono <- matrix(0,nrow=nG,ncol=nSim)
for (nsim in 1:nSim){
  cat(nsim,"out of",nSim,"\n")
  nSimMt <- round(rowSums(V)/dwn)
  for (i in 1:nG){
    set.seed(NULL)
    sim <- sample( x=1:nMtTps , size=nSimMt[i], replace=TRUE, prob=V[i,])
    sampleV[i,] <- tabulate(sim,nbins=nMtTps)
  }
  ResultFixTri <- NMFglmSQR(Data = sampleV, NoSignatures = 4, Signatures = TriRes$Signatures,tol=0.2,Seeds=sample(1:100,3))
  #ResultFixMix <- NMFglmSQR(Data = sampleV, NoSignatures = 4, Signatures = MixRes$Signatures,tol=0.2,Seeds=sample(1:100,3))
  ResultFixDi <- NMFglmSQR(Data = sampleV, NoSignatures = 4, Signatures = DiRes$Signatures,tol=0.2,Seeds=sample(1:100,3))
  ResultFixMono <- NMFglmSQR(Data = sampleV, NoSignatures = 4, Signatures = MonoRes$Signatures,tol=0.2,Seeds=sample(1:100,3))
  
  ## Compare exposures: Cosine similarity for each patient
  ResCosineMatTri[,nsim] = cosMatch(TriRes$Exposures,ResultFixTri$Exposures)$cossim
  ResCosineMatDi[,nsim] = cosMatch(DiRes$Exposures,ResultFixDi$Exposures)$cossim
  ResCosineMatMono[,nsim] = cosMatch(MonoRes$Exposures,ResultFixMono$Exposures)$cossim
}

##------------------------------------------------------------------
## Compare true exposures and estimated exposures
##------------------------------------------------------------------
## Mean
mnTri <- rowMeans(ResCosineMatTri)
mnMix <- rowMeans(ResCosineMatMix)
mnDi <- rowMeans(ResCosineMatDi)
mnMono <- rowMeans(ResCosineMatMono)
# plot mean
plot(mnTri,ylim=c(0.5,1),col="red",pch=19,cex=1)
points(mnMix,col="blue",pch=19,cex=0.8)
points(mnDi,col="green",pch=19,cex=0.6)
points(mnMono,col="orange",pch=19,cex=0.5)
## Var
vrTri <- apply(ResCosineMatTri,1,var)
vrMix <- apply(ResCosineMatMix,1,var)
vrDi <- apply(ResCosineMatDi,1,var)
vrMono <- apply(ResCosineMatMono,1,var)
plot(vrTri,col="red",pch=19,cex=0.8)
points(vrMix,col="blue",pch=19,cex=0.6)
points(vrMix,col="blue",pch=19,cex=0.6)
## Quantiles
qTri <- apply(ResCosineMatTri,1,quantile,probs=c(0.05,0.95))
qMix <- apply(ResCosineMatMix,1,quantile,probs=c(0.05,0.95))
## Mean with quantiles

plot(1:21,mnTri,ylim=c(0,1),col="red",pch=19,cex=0.9,
     main=paste("Ratio of down sampling:",dwn),
     xlab="Patient",ylab="Cosine Similarity",cex.lab=1.2)
axis(side=1, at=0:21, labels = TRUE)
points(1:21+0.3,mnMix,col="blue",pch=19,cex=0.9)
for (i in 1:21){
  points(rep(i,2),qTri[,i],col="red",type="l",lwd=2)
  points(rep(i+0.3,2),qMix[,i],col="blue",type="l",lwd=2)
}
lines(c(14,15),rep(0.5,2),col="red",lwd=2)
text(15,0.5,"Sole tri-nucleotide",pos=4,cex=1.2)
lines(c(14,15),rep(0.4,2),col="blue",lwd=2)
text(15,0.4,"Mixture",pos=4,lwd=2,cex=1.2)
points(14.5,0.3,pch=19,col="black",cex=1.2)
text(15,0.3,"Mean",pos=4,cex=1.2)
lines(c(14,15),rep(0.2,2),lty=1,lwd=2)
text(15,0.2,"Quantile interval",pos=4,cex=1.2)
text(1:21,rep(0.05,15),nSimMt,cex=0.8)
text(0.5,0.12,"Number of sampled mutations",pos=4,cex=1.2)

ResultFixDi <- NMFglmSQR(Data = V, NoSignatures = 4, Signatures = DiRes$Signatures,tol=0.2,Seeds=3)
cosMatch(DiRes$Signatures,ResultFixDi$Signatures)
