########################################################################
## Simulation study with iter and init to reach same minimum for Di 
#######################################################################
#rm(list=ls())
setwd("~/projects/paramNMF_ms")

source("BRCA/loadBRCA214models.R")

DiRes = resFactors[[11]]

##-----------------------------------------------------------------------
## simulation for same data
##-----------------------------------------------------------------------
sampleV <- V
nG <- dim(V)[1]     # Number of patients (genomes)
nMtTps <- dim(V)[2] # Number of mutation types
nMt <- rowSums(V)
noSig = 4

nSim <- 10  # 5 if small study. 100 if intermediate. 200 if large.
init = c(100,500,1000)   # number of initialisations
smallIter = c(10,50,100)
tol = 1e-8

resultsDi = matrix(list(), nrow = 3, ncol = 3)
rownames(resultsDi) = smallIter
colnames(resultsDi) = init

countResult = matrix(0, nrow = 3, ncol = 3)
rownames(countResult) = smallIter
colnames(countResult) = init
 for(x in smallIter){
   for(y in init){
     cat("Iterations :", x, "Initializations:", y, "\n")
     for(nsim in 1:nSim){
       cat(nsim,"out of",nSim,"\n")
     ResultFixDi <- nmfprm(data = sampleV, noSignatures = noSig, designMatrices = rep(list(Mdi),noSig),tol=tol,initial = y, smalliter = x)
     # Compare signatures: Cosine similarity for each patient
     resultsDi[paste0(x),paste0(y)][[1]][nsim] = ResultFixDi$gkl
     if(abs(ResultsFixDi$gkl - DiRes$gkl) < 2){
       countResult[paste0(x),paste0(y)] = countResult[paste0(x),paste0(y)] + 1
     }
     }
     save(resultDi,countResult, file = "DisimulationBRCA214.RData")
     cat(countResult,"\n")
   }
}
print(countResult)
