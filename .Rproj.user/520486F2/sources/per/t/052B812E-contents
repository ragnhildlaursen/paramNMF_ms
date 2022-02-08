setwd('U:/Projects/NMFmodelSelection/ModelSelectionGLM/ms_glm')
source("modelselection.R")
source('BRCA21MonoDiTriCompare.R')
TriExposures <- ResultTri$Exposures[c(3,4,2,1),]
MixExposures <- ResultMix$Exposures[c(1,3,2,4),]
##-----------------------------------------------------------------------
## Sample a number of mutations from each patient (downsampling)
##-----------------------------------------------------------------------
sampleV <- V
nMtTps <- dim(V)[1] # Number of mutation types
nG <- dim(V)[2]     # Number of patients (genomes)
nMt <- colSums(V)
dwn <- 50  # Down sampling 100, 50, 20, 10
nSim <- 100  # 5 if small study. 100 if intermediate. 200 if large.
ResCosineMatTri <- matrix(0,nrow=nG,ncol=nSim) 
ResCosineMatMix <- matrix(0,nrow=nG,ncol=nSim) 
for (nsim in 1:nSim){
  cat(nsim,"out of",nSim,"\n")
  nSimMt <- round(colSums(V)/dwn)
  for (i in 1:nG){
    set.seed(NULL)
    sim <- sample( x=1:nMtTps , size=nSimMt[i], replace=TRUE, prob=V[,i])
    sampleV[,i] <- tabulate(sim,nbins=nMtTps)
  }
  ResultFixTri <- NMFglmSQRfixedSignature(sampleV,TriSignatures,tol=0.2,Seeds=1)
  ResultFixMix <- NMFglmSQRfixedSignature(sampleV,MixSignatures,tol=0.2,Seeds=1)
  ## Compare exposures: Cosine similarity for each patient
  for (i in 1:nG){
    a <- TriExposures[,i]
    b <- ResultFixTri$Exposures[,i]
    ResCosineMatTri[i,nsim] <- sum(a*b)/sqrt(sum(a^2)*sum(b^2))
    a <- MixExposures[,i]
    b <- ResultFixMix$Exposures[,i]
    ResCosineMatMix[i,nsim] <- sum(a*b)/sqrt(sum(a^2)*sum(b^2))
  }
}
##------------------------------------------------------------------
## Compare true exposures and estimated exposures
##------------------------------------------------------------------
plot(ResCosineMatTri[,1],type="l",col="red",ylim=c(0,1))
points(ResCosineMatMix[,1],type="l",col="blue")
for (nsim in 1:nSim){
  points(ResCosineMatTri[,nsim],type="l",col="red",ylim=c(0,1))
  points(ResCosineMatMix[,nsim],type="l",col="blue")
}
## Mean
mnTri <- rowMeans(ResCosineMatTri)
plot(mnTri,ylim=c(0,1),col="red",pch=19,cex=0.8)
mnMix <- rowMeans(ResCosineMatMix)
points(mnMix,col="blue",pch=19,cex=0.6)
## Var
vrTri <- apply(ResCosineMatTri,1,var)
plot(vrTri,col="red",pch=19,cex=0.8)
vrMix <- apply(ResCosineMatMix,1,var)
points(vrMix,col="blue",pch=19,cex=0.6)
## Quantiles
qTri <- apply(ResCosineMatTri,1,quantile,probs=c(0.05,0.95))
qMix <- apply(ResCosineMatMix,1,quantile,probs=c(0.05,0.95))
## Mean with quantiles
#pdf(file=paste("DownSampling",dwn,".pdf",sep=""),width=8,height=6)
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
#dev.off()
##---------------------------------------------------------------------
## Estimate exposures for fixed signatures 
##---------------------------------------------------------------------
## SQUAREM version
library('SQUAREM')
NMFglmSQRfixedSignature = function(Data, Signatures, tolerance = 1e-2, 
                     maxIter = 5000, Seeds = c(1,2,3)){
  NoSignatures <- ncol(Signatures)
  MutationTypes = dim(Data)[1]  # mutation types
  Genomes       = dim(Data)[2]  # genomes
  
  GKLvalues = rep(0,length(Seeds)) # vector of different Generalised Kullback Leibler(GKL) values
  #Signaturelist = list()            # list of signature matrices
  Exposurelist = list()             # list of exposure matrices
  
  ##Function with one E and M step
  EMstep = function(x){
    par = exp(x)
    #Signatures = matrix(par[1:(MutationTypes*NoSignatures)], nrow = MutationTypes, ncol = NoSignatures)
    Exposures = matrix(par, nrow = NoSignatures, ncol = Genomes)
    
    #EstimateOfData = Signatures%*%Exposures
    
    #regularUpdate = as.vector(Signatures * ((Data/EstimateOfData) %*% t(Exposures)))
    #Signatures = matrix(glm.update(regularUpdate,DesignMatrix), ncol = NoSignatures) # glm update of Signatures
    
    EstimateOfData = Signatures%*%Exposures
    Exposures = Exposures * (t(Signatures) %*% (Data/EstimateOfData))# update of Exposure
    #Exposures = diag(1/rowSums(Exposures)) %*% Exposures             # make sure the column sum to one
    
    par = as.vector(Exposures)
    par[par <= 0] = 1e-10
    logpar = log(par)
    return(logpar)
  }
  
  # Function to create GKL value
  gklobj = function(x){
    par = exp(x)
    #Signatures = matrix(par[1:(MutationTypes*NoSignatures)], nrow = MutationTypes, ncol = NoSignatures)
    Exposures = matrix(par, nrow = NoSignatures, ncol = Genomes)
    
    EstimateOfData = Signatures%*%Exposures
    GKL <- gkl.dev(as.vector(Data),as.vector(EstimateOfData)) # GKLD value
    return(GKL)
  }
  
  for(i in 1:length(Seeds)){
    set.seed(Seeds[i])                                              # Setting fixed seed 
    #Signatures = matrix(runif(NoSignatures*MutationTypes), 
    #                    ncol = NoSignatures, nrow = MutationTypes)  # Initialize Signatures
    Exposures = matrix(runif(Genomes*NoSignatures), 
                       ncol = Genomes, nrow = NoSignatures)         # Initialize Exposures
    #Exposures = diag(1/rowSums(Exposures)) %*% Exposures            # make sure the rows sum to one
    
    #Initial = c(as.vector(Signatures),as.vector(Exposures))
    Initial = as.vector(Exposures)
    #SQUAREM run of the EM algorithm
    ResultSqr = squarem(Initial, fixptfn = EMstep, objfn = gklobj, control = list(tol = tolerance, maxiter = maxIter))
    par = exp(ResultSqr$par) # parameters
    #Signatures = matrix(par[1:(MutationTypes*NoSignatures)], nrow = MutationTypes, ncol = NoSignatures)
    Exposures = matrix(par, nrow = NoSignatures, ncol = Genomes)
    
    GKLvalues[i] = gklobj(log(par))
    #Signaturelist[[i]] = Signatures
    Exposurelist[[i]] = Exposures
  }
  optimal = which.min(GKLvalues)
  #SignatureOptimal = Signaturelist[[optimal]]
  ExposureOptimal = Exposurelist[[optimal]]
  
  # Make columns of signatures sum to one
  #ExposureOptimal = diag(colSums(SignatureOptimal)) %*% ExposureOptimal
  #SignatureOptimal = SignatureOptimal %*% diag(1/colSums(SignatureOptimal))
  
  Output = list()
  #Output$Signatures = SignatureOptimal
  Output$Exposures = ExposureOptimal
  Output$gkl = GKLvalues[optimal]
  
  return(Output)
}

