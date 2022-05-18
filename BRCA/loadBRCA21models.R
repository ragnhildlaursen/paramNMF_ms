#######################################
## Load BRCA21 results and models
######################################
setwd("~/projects/paramNMF_ms/")

source("ModelSelection.R")
library(Rcpp)
library(RcppArmadillo)
sourceCpp("NMF2.cpp")

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
Mfull = model.matrix(~L*M*R)      # full model
Mdi = model.matrix(~L*M + M*R)    # di-nucleotide model
Mmono = model.matrix(~L + M + R)  # multiplicative model

load("BRCA/result/BRCA21modelFactors4sig500initv2.RData")
load("BRCA/result/BRCA21modelsummary4sig500initv2.RData")

