setwd('U:/Projects/NMFmodelSelection/ModelSelectionGLM/ms_glm')
## Include core functions
source("modelselection.R")
## Load the data
load("BRCA21.RData")
## The different factors 
left = c(rep("lA",4), rep("lC",4), rep("lG",4), rep("lT",4))
right = c(rep(c("rA","rC","rG","rT"),4))
mut = c("C>A","C>G","C>T","T>A","T>C","T>G")
## The factors for 96 different mutation types 
R = factor(rep(right,6))
L = factor(rep(left,6))
M = factor(rep(mut, each = 16))
## Model matrices
M0 = model.matrix(~0+L*M*R)         # tri-nucleotide signature
M1 = model.matrix(~0+L*M + M*R)     # di-nucleotide signature
M2 = model.matrix(~0+L + M + R)     # mono-nucleotide signature
##-----------------------------------------
## Estimate the three models:
## K=4 mono-nucleotide
## K=4 di-nucleotide
## K=4 tri-nucleotide
##-----------------------------------------
## Mono-nucleotide model 
library("magic")
DesignMatrix = adiag(M2,M2,M2,M2) 
ResultMono = NMFglmSQR(Data=V, NoSignatures=4,DesignMatrix = DesignMatrix, 
                    tolerance = 0.2, Seeds = 20) # Seeds=20 gives minimum!
cat(ResultMono$gkl) # 4647.098
## Tri-nucleotide model
DesignMatrix = adiag(M0,M0,M0,M0) 
ResultTri = NMFglmSQR(Data=V, NoSignatures=4,DesignMatrix = DesignMatrix, 
                    tolerance = 0.1, Seeds = 18) # Seeds=18 gives minimum!
cat(18,ResultTri$gkl,"\n") # 1474.534
## Di-nucleotide model 
DesignMatrix = adiag(M1,M1,M1,M1) 
ResultDi = NMFglmSQR(Data=V, NoSignatures=4,DesignMatrix = DesignMatrix, 
                    tolerance = 0.1, Seeds = 24) # Seeds=24 gives minumum!
cat(ResultDi$gkl,"\n") # 2070.915
## Mixture model
DesignMatrix = adiag(M0,M1,M1,M2)
ResultMix = NMFglmSQR(Data=V, NoSignatures=4,DesignMatrix = DesignMatrix, 
                     tolerance = 0.1, Seeds = 13) # Seeds=13 gives minumum!
cat(ResultMix$gkl,"\n") # 1860.876
#for (i in 40:50){
#  DesignMatrix = adiag(M1,M1,M1,M1)
#  Result = NMFglmSQR(Data=V, NoSignatures=4,DesignMatrix = DesignMatrix, 
#                       tolerance = 0.2, Seeds = i) 
#  cat(Result$gkl,"\n") 
#}
##----------------------------------------------------------
## Define cosine similarity
## and note their nice similarity
CosSim <- function(a,b) {sum(a*b)/sqrt(sum(a^2)*sum(b^2))}
##----------------------------------------------------------------------
## Align the signatures for the three models according to 
## Figure 4 in Alexandrov et al (2013, Cell) 
## (Camilla Kudahl master thesis page 46, Appendix A)
##-----------------------------------------------------------------------
TriSignatures <- ResultTri$Signatures[,c(3,4,2,1)]
DiSignatures <- ResultDi$Signatures[,c(1,3,2,4)]
MonoSignatures <- ResultMono$Signatures[,c(4,1,2,3)]
MixSignatures <- ResultMix$Signatures[,c(1,3,2,4)]
MixTypes <- c("(tri-nucleotide)","(di-nucleotide)","(di-nucleotide)","(mono-nucleotide)")
for (i in 1:4){
  cat("Signature",i,"\n")
  cat("Mono versus Di",CosSim(MonoSignatures[,i],DiSignatures[,i]),"\n")
  cat("Mono versus Tri",CosSim(MonoSignatures[,i],TriSignatures[,i]),"\n")
  cat("Di versus Tri",CosSim(DiSignatures[,i],TriSignatures[,i]),"\n")
  cat("Mix versus Mono",CosSim(MixSignatures[,i],MonoSignatures[,i]),"\n")
  cat("Mix versus Di",CosSim(MixSignatures[,i],DiSignatures[,i]),"\n")
  cat("Mix versus Tri",CosSim(MixSignatures[,i],TriSignatures[,i]),"\n")
}
i <- 1
ymx <- c(0.2,0.5,0.2,0.5)
plot(1:96,MonoSignatures[,i],ylim=c(0,ymx[i]),type="h",col="orange",
     xlab="Mutation type",ylab="Probability")
points(1:96+0.2,DiSignatures[,i],type="h",col="blue")
points(1:96+0.4,TriSignatures[,i],type="h",col="red")
points(1:96+0.6,MixSignatures[,i],type="h",col="green")
## Plot and compare the signatures
specifyDecimal <- function(x, k) trimws(format(round(x, k), nsmall=k))
for (i in 1:4){
#  pdf(file=paste("CmprSgntr",i,".pdf",sep=""),width=10,height=5)
  ymx <- c(0.2,0.5,0.2,0.5)
  plot(1:96,MonoSignatures[,i],ylim=c(0,ymx[i]),type="h",col="orange",
       xlab=ifelse(i==4,"Mutation type",""),ylab="Probability",xaxt="n")
  if (i==4){
    axis(side=1, at=c(1,16,32,48,64,80,96),labels=TRUE)
    axis(side=1, at=c(8,24,40,56,72,88),labels=mut,tick=FALSE)
  }
  points(1:96+0.2,DiSignatures[,i],type="h",col="darkblue")
  points(1:96+0.4,TriSignatures[,i],type="h",col="red")
  points(1:96+0.6,MixSignatures[,i],type="h",col="darkgreen")
  legend("topright",c("Exclusive mono-nucleotide",
                      "Exclusive di-nucleotide",
                      "Exclusive tri-nucleotide",
                      paste("Mixture",MixTypes[i])),
         lty=1,col=c("orange","darkblue","red","darkgreen"),bty="n",cex=1.2,lwd=2)
  text(0,ymx[i]*0.97,paste("Signature",i),pos=4,cex=1.5)
  ## Cosine similarity
  text(0,ymx[i]*0.89,"Cosine Similarities:",pos=4,cex=1.2)
  text(0,ymx[i]*0.82,"Mono and Di",pos=4) 
  text(17,ymx[i]*0.82,specifyDecimal(CosSim(MonoSignatures[,i],DiSignatures[,i]),2),pos=4) 
  text(0,ymx[i]*0.77,"Mono and Tri",pos=4) 
  text(17,ymx[i]*0.77,specifyDecimal(CosSim(MonoSignatures[,i],TriSignatures[,i]),2),pos=4)
  text(0,ymx[i]*0.72,"Di and Tri",pos=4) 
  text(17,ymx[i]*0.72,specifyDecimal(CosSim(DiSignatures[,i],TriSignatures[,i]),2),pos=4) 
  text(0,ymx[i]*0.67,"Mix and Mono",pos=4) 
  text(17,ymx[i]*0.67,specifyDecimal(CosSim(MixSignatures[,i],MonoSignatures[,i]),2),pos=4) 
  text(0,ymx[i]*0.62,"Mix and Di",pos=4) 
  text(17,ymx[i]*0.62,specifyDecimal(CosSim(MixSignatures[,i],DiSignatures[,i]),2),pos=4)
  text(0,ymx[i]*0.57,"Mix and Tri",pos=4) 
  text(17,ymx[i]*0.57,specifyDecimal(CosSim(MixSignatures[,i],TriSignatures[,i]),2),pos=4) 
#  dev.off()
}
