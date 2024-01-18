# paramNMF_ms
This project parametrizes mutational signatures to give more robust and interpretable results. 

The data and codes here can be used to recreate the results seen in the paper 'Robust estimation of mutational signatures by adaptive learning of nucleotide interaction terms in non–negative matrix factorization' by Ragnhild Laursen, Lasse Marretty and Asger Hobolth. 

The code files *ModelSelection.R* and *NMF2.cpp* are the essential code used to run non-negative matrix factorization with parametrization of the signatures. The ModelSelection.R file include the original codes and NMF2.cpp includes a main function **nmfprm()** to find the estimates, where it is implemented i C++ to achieve faster convergence. Are you not able to use C++, then you can use the original function in *ModelSelection.R* called **NMFglmSQR()**.

## Analysis of BRCA datasets 
The code for recreating the results for the BRCA data sets can be found in the **BRCA** folder. This folder consist of both data sets, codes to create the results and a folder that holds the results. 

The two breast cancer datasets of mutational counts are *BRCA21.RData* and *BRCA214.RData*. They consist of a datamatrix of size 21 x 96 and 214 x 96, respectively.

The coding files do the folowing:
 - *loadBRCA21models.R* - loads the BRCA21 dataset and creates the factors and model matrices for the parametrization of the signatures.
 - *loadBRCA214models.R* - loads the BRCA214 dataset and creates the factors and model matrices for the parametrization of the signatures.
 - *optimalSignaturesBRCA.R* - recovers the fit and estimates of $W$ and $H$ for each of the parametrized models. 
 - *RobustExposures.R* - find the cosine similarity between the orignal exposures and the reestimated ones from the downsampled data. 
 - *RobustSignatures.R* - find the cosine similarity between the original signatures and the reestimed ones after applying parametric bootstrapping to the original results.

## Analysis of UCUT datasets

The folder **UCUT** contains all the necessary data and code to recreate the results for the UCUT analysis. The structure of the folder is equivalent to the **BRCA** folder, where the dataset is the file *UCUT_5_all.RData*. The data matrix here is of size 26 x 1536.


## Code for plotting the figures 
The code for plotting the figures in the paper are found in *plotcode* and the resulting plots can be found in the folder *plots*.

## Example for using the function and model

```r
# set working directory
setwd("~/paramNMF_ms/")

# libraries used
library(Rcpp)
library(RcppArmadillo)

# loading the functions
source("ModelSelection.R")
sourceCpp("NMF2.cpp")

# load count data, where the rows are different patients and the columns correspond to the different mutation types
load("BRCA/BRCA214.RData")

##--------------------------------------------------------
## Factor variables
##--------------------------------------------------------
## The factors for the different mutation types 
L = factor(substr(colnames(V), start = 1, stop = 1)) # a vector of the left flanking nucleotide
M = factor(substr(colnames(V), start = 3, stop = 5)) # a vector of the base mutation
R = factor(substr(colnames(V), start = 7, stop = 7)) # a vector of the right flanking nucleotide

##--------------------------------------------------------
## Parametrizations of a signature
##--------------------------------------------------------
## Model matrices
Mfull = model.matrix(~L*M*R)      # full model
Mdi = model.matrix(~L*M + M*R)    # di-nucleotide model
Mmono = model.matrix(~L + M + R)  # mono-nucleotide model

##--------------------------------------------------------
## Run a mixture model, where the first signatures follows mono-nucleotide parametrization, second a di-nucleotide parametrization and the third with the full/regular parametrization.
##--------------------------------------------------------

res1 <- nmfprm(data=V,noSignatures=3,designMatrices=list(Mmono,Mdi,Mdi,Mfull))

# The exposures for the different signatures
W = res1$exposures

# The signatures that will follow their respective parametrizations defined in the function eg. the first signature will follow the mono-nucleotide parametrization.
H = res1$signatures

##--------------------------------------------------------
## Run a mode, where all the signatures follow a di-nucleotide parametrization
##--------------------------------------------------------

res2 <- nmfprm(data=V,noSignatures=8,designMatrices=rep(list(Mdi),8))

```

