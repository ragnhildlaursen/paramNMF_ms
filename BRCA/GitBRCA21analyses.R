##------------------------------------------------------------
## BRCA21 data analysis with number of signatures equal to four
##------------------------------------------------------------
## We consider 15 models with the three signatures 
## Mono-nucleotide (Shiraishi et al), 
## Di-nucleotide (our proposed signature), 
## Tri-nucleotide (Alexandrov et al)
## Run each model 30 times with different initializations 
##------------------------------------------------------------
source("GitModelSelection.R")

# load BRCA data 
load("BRCA/BRCA21.RData")
V = t(V) # rows should be patients and columns mutation types
##--------------------------------------------------------
## Factors
##--------------------------------------------------------
## The factors for 96 different mutation types 
L = factor(substr(colnames(V), start = 1, stop = 1))
M = factor(substr(colnames(V), start = 3, stop = 5))
R = factor(substr(colnames(V), start = 7, stop = 7))

##--------------------------------------------------------
## Parametrizations of a signature
##--------------------------------------------------------
## Model matrices
Mfull = model.matrix(~0+L*M*R)      # full model
Mdi = model.matrix(~0+L*M + M*R)    # di-nucleotide model
Mmono = model.matrix(~0+L + M + R)  # multiplicative model
models = list(Mmono,Mdi,Mfull)      # list of 3 models

## List of 15 models for the 3 parametrizations with 4 signatures
noSignatures=4
ModelCombinations = combinations(length(models), noSignatures, 
                                 repeats.allowed = TRUE)
nModels <- nrow(ModelCombinations)
MList = lapply(1:nModels,function(x) lapply(ModelCombinations[x,], function(y) models[[y]]))

##-----------------------------------------------------------
## Matrix with summary of the results
##-----------------------------------------------------------
resMat <- matrix(0,nrow=nModels,ncol=7)
colnames(resMat) <- c("nprm1","nprm2","nprm3","nprm4","nprmtot","seed","GKL")
resFactors = list()
## First run the EM algorithm 30 times with high tolerance to 
## identify global minimum 
## Second run EM algorithm with seed from global minimum but 
## with low tolerance (i.e. until convergence)
## Here we use high tolerance=25 and low tolerance=5
## to run the algorithm in a few minutes.
## In order to obtain the results in the paper you *must* run with
## high tolerance=1 and low tolerance=0.2
## On a normal latop this takes around 5 minutes. 
## The GKL is summarized in resMat. 
## We provide resMat for our run in file "GitUCUT26results.txt"
high.tolerance <- 0.5
low.tolerance <- 0.01
repshigh = 30
for (m in 1:nModels){
  nprm1 <- ncol( MList[[m]][[1]] )
  nprm2 <- ncol( MList[[m]][[2]] )
  nprm3 <- ncol( MList[[m]][[3]] )
  nprm4 <- ncol( MList[[m]][[4]] )
  resMat[m,1] <- nprm1
  resMat[m,2] <- nprm2
  resMat[m,3] <- nprm3
  resMat[m,4] <- nprm4
  resMat[m,5] <- nprm1+nprm2+nprm3+nprm4
  ## List of models chosen
  ## Run the EM algorithm many times with a high tolerance to identify global minimum
  tmp <- rep(0,repshigh)
  seed = sample(1:10000,repshigh)
  cat("Model:",m,"","\n")
  cat("EM Run: ")
  for (j in 1:length(tmp)){
    cat(j,"")
    tmp[j] = NMFglmSQR(Data=V, NoSignatures=4,DesignMatrix = MList[[m]], 
                       tolerance =high.tolerance, Seeds = seed[j])$gkl
  }
  ## Choose the best initial value and run again until convergence, i.e.
  ## stop when the tolerance is small
  cat("\n")
  jmin <- which.min(tmp)
  resMat[m,6] <- jmin
  res <- NMFglmSQR(Data=V,NoSignatures=4,DesignMatrix=MList[[m]],
                   tolerance=low.tolerance,Seeds=seed[jmin])
  resFactors[[m]] = res
  
  ## print results
  resMat[m,7] <- res$gkl
  cat("Final result:","nprm:",nprm1,nprm2,nprm3,nprm4,
      "; seed:",seed[jmin],", GKL:",res$gkl,"\n" )
}
#resMat <- resMat[sort(resMat[,"nprmtot"],index=TRUE)$ix,]
print(resMat)
par(mfrow = c(1,1))
BIC = 2*resMat[,"GKL"] + log(21*96)*resMat[,"nprmtot"]
plot(BIC)
which.min(BIC)
plot(resMat[,"GKL"])

save(resFactors, file = "BRCAmodelFactors.RData")
save(resMat, file = "BRCAmodelRes.RData")
