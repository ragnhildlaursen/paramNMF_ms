#######################################
## Load BRCA21 results and models
######################################
setwd("~/projects/paramNMF_ms/")

source("GitModelSelection.R")
library(Rcpp)
library(RcppArmadillo)
sourceCpp("fastercode/NMF2.cpp")

# load BRCA data 
load("BRCA/BRCA214.RData")
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

load("BRCA/result/BRCA214modelFactors4sig500init.RData")
load("BRCA/result/BRCA214modelsummary4sig500init.RData")

