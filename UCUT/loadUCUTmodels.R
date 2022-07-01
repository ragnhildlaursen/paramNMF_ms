#######################################
## Load BRCA results and models
######################################
library(Rcpp)
library(RcppArmadillo)
setwd("~/projects/paramNMF_ms/")
source("ModelSelection.R")
sourceCpp("NMF2.cpp")

## Load UCUT data
load("UCUT/UCUT_5_all.RData")

V = t(Vall)   # data (No. patients) x (No. Mutation types)
## UCUT has 1536 mutation types and 26 patients 
V[V == 0] = .Machine$double.eps # Zero entries replaced with small epsilon to avoid division by zero in EM-algorithm
##--------------------------------------------------------
## Factors
##--------------------------------------------------------
## First and second left flanking nucleotide 
l1 = factor(substr(colnames(V), start = 1, stop = 1))
l2 = factor(substr(colnames(V), start = 9, stop = 9))
## Actual point mutation
m = factor(substr(colnames(V), start = 3,stop = 5))
## First and second right flanking nucleotide 
r1 = factor(substr(colnames(V), start = 7, stop = 7))
r2 = factor(substr(colnames(V), start = 11, stop = 11))

##--------------------------------------------------------
## Parametrizations of a signature
##--------------------------------------------------------
Mmono = model.matrix(~0+l2+l1+m+r1+r2)          # Mono-nucleotide model
Mdi = model.matrix(~0+m*l2+l1*m+m*r1+m*r2)      # Di-nucleotide interaction with mutation
Mblend = model.matrix(~0+l2+l1*m+m*r1+r2)       # Blended mono-di-nucleotide model
Mcombi = model.matrix(~0+l2+l1*m*r1+r2)         # Combined mono-tri-nucleotide model
Mtri =  model.matrix(~0+l1*m*r1)                # Tri-nucleotide model (only one flanking nucleotide)
Mnghbr = model.matrix(~0+l2*l1+l1*m+m*r1+r1*r2) # Di-nucleotide interaction with neighbour
Mfull = model.matrix(~0+l2*l1*m*r1*r2)          # Full parametrized model

Mditri = model.matrix(~0+l2*l1 + l1*m*r1 + r1*r2)

## List of the 21 models
# MList <- list(list(Mmono,Mmono),
#               list(Mdi,Mdi),list(Mmono,Mdi),
#               list(Mblend,Mblend),list(Mmono,Mblend),list(Mdi,Mblend),
#               list(Mcombi,Mcombi),list(Mmono,Mcombi),list(Mdi,Mcombi),list(Mblend,Mcombi),
#               list(Mtri,Mtri),list(Mmono,Mtri),list(Mdi,Mtri),list(Mblend,Mtri),
#               list(Mcombi,Mtri),
#               list(Mnghbr,Mnghbr),list(Mmono,Mnghbr),list(Mdi,Mnghbr),list(Mblend,Mnghbr),
#               list(Mcombi,Mnghbr),list(Mtri,Mnghbr))
# nModels <- length(MList)


load("UCUT/result/UCUTmodelFactors.RData")
load("UCUT/result/UCUTmodelsummary.RData")

