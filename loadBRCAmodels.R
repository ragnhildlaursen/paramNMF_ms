#######################################
## Load BRCA results and models
######################################
source("GitModelSelection.R")

# load BRCA data 
load("BRCA/BRCA21.RData")
V = t(V) # change to dimension patients x mutation types
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

load("BRCAmodelFactors.RData")
TriRes = resFactors[[15]]
DiRes = resFactors[[11]]
MixRes = resFactors[[9]]
MonoRes = resFactors[[1]]