################################################
## Checking robustness of mutational signatures 
################################################
source("GitModelSelection.R")

# load BRCA data 
load("BRCA/BRCA21.RData")

##--------------------------------------------------------
## Factors
##--------------------------------------------------------
## The factors for 96 different mutation types 
L = factor(substr(rownames(V), start = 1, stop = 1))
M = factor(substr(rownames(V), start = 3, stop = 5))
R = factor(substr(rownames(V), start = 7, stop = 7))

##--------------------------------------------------------
## Parametrizations of a signature
##--------------------------------------------------------
## Model matrices
Mfull = model.matrix(~0+L*M*R)      # full model
Mdi = model.matrix(~0+L*M + M*R)    # di-nucleotide model
Mmono = model.matrix(~0+L + M + R)  # multiplicative model

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