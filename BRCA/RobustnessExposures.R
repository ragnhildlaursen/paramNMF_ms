rm(list=ls())
################################################################
## Robustness of exposures
################################################################
setwd("~/projects/paramNMF_ms")
source("BRCA/loadBRCAmodels.R")
glm.fit()
quasipoisson()
##-----------------------------------------------------------------------
## Sample a number of mutations from each patient (downsampling)
##-----------------------------------------------------------------------
sampleV <- V
nG <- dim(V)[1]     # Number of patients (genomes)
noSig = 4
nMtTps <- dim(V)[2] # Number of mutation types
nMt <- rowSums(V)
dwn <- 10  # Down sampling 100, 20, 10
nSim <- 100  # 5 if small study. 100 if intermediate. 200 if large.
tol = 0.1
ResCosineMatTri <- matrix(0,nrow=noSig,ncol=nSim) 
#ResCosineMatMix <- matrix(0,nrow=noSig,ncol=nSim) 
ResCosineMatDi <- matrix(0,nrow=noSig,ncol=nSim) 
ResCosineMatMono <- matrix(0,nrow=noSig,ncol=nSim)
for (nsim in 1:nSim){
  cat(nsim,"out of",nSim,"\n")
  nSimMt <- round(rowSums(V)/dwn)
  for (i in 1:nG){
    set.seed(NULL)
    sim <- sample( x=1:nMtTps , size=nSimMt[i], replace=TRUE, prob=V[i,])
    sampleV[i,] <- tabulate(sim,nbins=nMtTps)
  }
  ResultFixTri <- NMFglmSQR(Data = sampleV, NoSignatures = noSig, Signatures = TriRes$Signatures,tol=tol)
  #ResultFixMix <- NMFglmSQR(Data = sampleV, NoSignatures = 4, Signatures = MixRes$Signatures,tol=tol)
  ResultFixDi <- NMFglmSQR(Data = sampleV, NoSignatures = noSig, Signatures = DiRes$Signatures,tol=tol)
  ResultFixMono <- NMFglmSQR(Data = sampleV, NoSignatures = noSig, Signatures = MonoRes$Signatures,tol=tol)
  
  ## Compare exposures: Cosine similarity for each patient
  ResCosineMatTri[,nsim] = cosMatch(t(TriRes$Exposures),t(ResultFixTri$Exposures))$cossim
  ResCosineMatDi[,nsim] = cosMatch(t(DiRes$Exposures),t(ResultFixDi$Exposures))$cossim
  #ResCosineMatMix[,nsim] = cosMatch(t(MixRes$Exposures),t(ResultFixMix$Exposures))$cossim
  ResCosineMatMono[,nsim] = cosMatch(t(MonoRes$Exposures),t(ResultFixMono$Exposures))$cossim
}
save(ResCosineMatTri,ResCosineMatDi,ResCosineMatMono, file = "ExposureBRCAdwn10.RData")
##------------------------------------------------------------------
## Compare true exposures and estimated exposures
##------------------------------------------------------------------
datTri = data.frame(s = t(ResCosineMatTri), type = "tri")
#m = cosMatch(TriRes$Signatures, MixRes$Signatures)$match
#datMix = data.frame(s = t(ResCosineMatMix[m,]), type = "mix")
m = cosMatch(TriRes$Signatures, DiRes$Signatures)$match
datDi = data.frame(s = t(ResCosineMatDi[m,]), type = "di")
m = cosMatch(TriRes$Signatures, MonoRes$Signatures)$match
datMono = data.frame(s = t(ResCosineMatMono[m,]), type = "mono")

dat1 = rbind(datTri,datDi,datMono)
dat2 = reshape(dat1, varying = colnames(dat1)[1:4], direction = "long")
dat2$type = factor(dat2$type, levels = c("mono","di","mix","tri"))
dat2$time = factor(dat2$time)
dat2$downsample = 1/dwn*100

# remove mix
#dat2 = dat2[dat2$type != "mix",]
dat20 = dat2

datall = rbind(dat100,dat50,dat20)
##------------------------------------------------------------------
## Compare true exposures and estimated exposures
##------------------------------------------------------------------
par(mfrow = c(2,2))
p1 = ggplot(dat2[dat2$time == 1,], aes(x=type, y=s, fill=type)) +
  geom_violin(position=position_dodge(1)) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.9) +
  ggtitle("Signature 1")
p2 = ggplot(dat2[dat2$time == 2,], aes(x=type, y=s, fill=type)) +
  geom_violin(position=position_dodge(1)) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.9) +
  ggtitle("Signature 2")
p3 = ggplot(dat2[dat2$time == 3,], aes(x=type, y=s, fill=type)) +
  geom_violin(position=position_dodge(1)) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.9) +
  ggtitle("Signature 3")
p4 = ggplot(dat2[dat2$time == 4,], aes(x=type, y=s, fill=type)) +
  geom_violin(position=position_dodge(1)) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.9) +
  ggtitle("Signature 4")
ggarrange(p1,p2,p3,p4, nrow = 2, ncol = 2, common.legend = TRUE)

ggplot(datall, aes(x=type, y=s, fill=type)) +
  geom_violin(position=position_dodge(1), scale = "width") +
  facet_grid(rows = vars(time), cols = vars(downsample), scales = "free") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.3) +

  ylab("Signatures")
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
