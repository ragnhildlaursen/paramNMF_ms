##########################################################
## Model selection: Estimate parameters for 
## Poisson-Parametrized-NMF 
##########################################################

## Include package "SQUAREM"
library("gtools") ## Required for function 'combinations'
library('SQUAREM')
###########################################################
#' Generalized Kullback Leibler divergence
#'
#' @param y data 
#' @param mu estimate of data
#'
#' @return generalized kullback leibler divergence
#' @export
#'
gkl.dev = function(y, mu){
  if(any(is.na(y)) | any(is.na(mu))) stop("Input include NaN values.")
  
  r = mu
  p = which(y > 0)
  r[p] = (y * (log(y)- log(mu)) - y + mu)[p]
  return(sum(r))
}
###############################################################
#' Poisson GLM update 
#'
#' @param y data 
#' @param X design matrix
#' @param maxiter maximum number of iterations
#' @param epsilon change of deviance for stopping 
#'
#' @return estimate of data with GLM dependency
#' @export
#'
glm.update = function(y,X, maxiter = 40, epsilon = 1e-8){
  if(length(y) <= ncol(X)){
    return(y)
  }else{
    mu = y + runif(length(y),0,0.1)            # starting value
    dev.old = gkl.dev(y,mu)                    # start deviance
    for(d in 1:maxiter){
      z = log(mu) + (y-mu)/mu
      w = sqrt(mu)
      fit = .Call(stats:::C_Cdqrls, X*w, z*w, epsilon, FALSE) # solving least squares
      coef = fit$coefficients                                 # coefficients
      eta = as.vector(X%*%coef)
      mu = exp(eta)
      dev.new = gkl.dev(y,mu)
      if(2*abs(dev.old - dev.new)/(0.1 + abs(2*dev.new)) < epsilon){break}
      dev.old = dev.new
    } 
    return(mu)
  }
}
#####################################################
#' Parametric model updates to factorize a non-negative data matrix made faster with SQUAREM
#'
#' @param Data Data matrix i.e. genomes in rows and mutation types in columns
#' @param NoSignatures Rank of factorization i.e. number of signatures. Default is length of DesignMatrix
#' @param DesignMatrix List of design matrices to describe parametrization of the signatures. Default is identity matrices.
#' @param tolerance Stopping tolerance of the change of GKL
#' @param maxIter Maximum number of iterations
#' @param Seeds The seeds for the initialization
#'
#' @return A list of the matrices derived by the factorization and the corresponding generalized Kullback-Leibler
#' \itemize{
#'  \item Signature  - Non-negative matrix of dimension  NoSignatures x ncol(Data), with rows summing to one
#'  \item Exposure   - Non-negative matrix of dimension nrow(Data) x NoSignatures
#'  \item gkl - Smallest Value of the Generalized Kullback-Leibler
#'  }
#'  
#' @export
#'
NMFglmSQR = function(Data, NoSignatures = length(DesignMatrix), 
                     DesignMatrix = rep(list(diag(ncol(Data))),NoSignatures), 
                     tolerance = 1e-2, maxIter = 5000, Seeds = c(1,2,3), Exposures = NULL, Signatures = NULL){
  
  if(!is.list(DesignMatrix)) stop("DesignMatrix needs to be a list of matrices.")
  if(NoSignatures != length(DesignMatrix)) stop("NoSignatures different from number of specified matrices in DesignMatrix.")
  if(any(sapply(DesignMatrix,nrow) != ncol(Data))) stop("The number of rows for matrices in DesignMatrix needs to equal the columns in Data.")
  
  fixSig = FALSE
  fixExp = FALSE
  if(!is.null(Signatures)){
    fixSig = TRUE
    if(any(dim(Signatures) != c(NoSignatures,ncol(Data)))){
      stop("The dimensions of the Signatures does not match the Data and NoSignatures") 
    }
  }
  if(!is.null(Exposures)){
    fixExp = TRUE
    if(any(dim(Exposures) != c(nrow(Data),NoSignatures))){
      stop("The dimensions of the Exposures does not match the Data and NoSignatures") 
  }
  }
  
  Genomes       = dim(Data)[1]  # genomes
  MutationTypes = dim(Data)[2]  # mutation types
  
  GKLvalues = rep(0,length(Seeds)) # Vector of different Generalised Kullback Leibler(GKL) values
  Signaturelist = list()           # list of signature matrices
  Exposurelist = list()            # list of exposure matrices
  
  ## Function with one E and M step
  EMstep = function(x, param = FALSE){
    par = exp(x)
    if(!fixExp){
      Exposures = matrix(par[c(1:(Genomes*NoSignatures))], nrow = Genomes, ncol = NoSignatures)
    }
    if(!fixSig){
      Signatures = matrix(par[-c(1:(Genomes*NoSignatures))], nrow = NoSignatures, ncol = MutationTypes)
    }
    
    if(!fixSig){
      EstimateOfData = Exposures%*%Signatures
      regularUpdate = Signatures * (t(Exposures) %*% (Data/EstimateOfData))
      if(param){
        glmUpdate = t(sapply(1:NoSignatures, function(k) glm.update(regularUpdate[k,], DesignMatrix[[k]]))) # glm update of Signatures
        Signatures = glmUpdate
      }else{
        Signatures = regularUpdate
      }
      Signatures = diag(1/rowSums(Signatures))%*%Signatures          # make sure the rows sum to one
    }
    
    if(!fixExp){
      EstimateOfData = Exposures%*%Signatures
      Exposures = Exposures * ((Data/EstimateOfData) %*% t(Signatures)) # update of exposures
    }
    
    
    par = c(as.vector(Exposures),as.vector(Signatures))
    par[par <= 0] = .Machine$double.eps
    logpar = log(par) # using log-scale to also allow negative values
    return(logpar)
  }
  
  # Function to create GKL value
  gklobj = function(x){
    par = exp(x)
    if(!fixExp){
      Exposures = matrix(par[c(1:(Genomes*NoSignatures))], nrow = Genomes, ncol = NoSignatures)
    }
    if(!fixSig){
      Signatures = matrix(par[-c(1:(Genomes*NoSignatures))], nrow = NoSignatures, ncol = MutationTypes)
    }
    
    EstimateOfData = Exposures%*%Signatures
    GKL <- gkl.dev(as.vector(Data),as.vector(EstimateOfData)) # GKLD value
    return(GKL)
  }
  
  for(i in 1:length(Seeds)){
    set.seed(Seeds[i])                                              # Setting fixed seed 
    if(!fixExp){
      Exposures = matrix(runif(Genomes*NoSignatures), 
                         nrow = Genomes, ncol = NoSignatures)       # Initialize Exposures
    }
    if(!fixSig){
    Signatures = matrix(runif(NoSignatures*MutationTypes), 
                        nrow = NoSignatures, ncol = MutationTypes)  # Initialize Signatures
    }
    
    Initial = c(as.vector(Exposures),as.vector(Signatures))
    
    #SQUAREM run of the EM algorithm
    ResultSqrFull = squarem(par = Initial, fixptfn = EMstep, objfn = gklobj, control = list(tol = tolerance, maxiter = 100))
    ResultSqr = squarem(par = ResultSqrFull$par, fixptfn = function(x) EMstep(x,param = T), objfn = gklobj, control = list(tol = tolerance, maxiter = maxIter))
    print(ResultSqr$fpevals)
    par = exp(ResultSqr$par) # parameters
    
    if(!fixExp){
      Exposures = matrix(par[c(1:(Genomes*NoSignatures))], nrow = Genomes, ncol = NoSignatures)
    }
    if(!fixSig){
      Signatures = matrix(par[-c(1:(Genomes*NoSignatures))], nrow = NoSignatures, ncol = MutationTypes)
    }
    GKLvalues[i] = gklobj(log(par))
    Signaturelist[[i]] = Signatures
    Exposurelist[[i]] = Exposures
  }
  optimal = which.min(GKLvalues)
  SignatureOptimal = Signaturelist[[optimal]]
  ExposureOptimal = Exposurelist[[optimal]]
  
  # Make columns of signatures sum to one
  ExposureOptimal =  ExposureOptimal%*%diag(rowSums(SignatureOptimal))
  SignatureOptimal = diag(1/rowSums(SignatureOptimal))%*%SignatureOptimal 
  
  Output = list()
  Output$Signatures = SignatureOptimal
  Output$Exposures = ExposureOptimal
  Output$gkl = GKLvalues[optimal]
  Output$prmSignatures = sapply(DesignMatrix,ncol)
  
  return(Output)
}

##################################################
## Cosine similarity between a vector x and y
###############################################
similarity <- function(x,y){
  dot.prod <- sum(x*y) 
  norm.x <- sqrt(sum(x^2))
  norm.y <- sqrt(sum(y^2))
  frac <- dot.prod / (norm.x * norm.y)
  return(as.numeric(frac))
}

#' @title Cosine similarity of mutational signatures
#'
#'
#' @param H1 Numeric matrix of mutational signatures. Each row should represent a signature.
#' @param H2 Numeric matrix of mutational signatures. Each row should represent a signature.
#'
#'
#' @return List of match and cosine similarities
#' - \texttt{match} Vector of indexes to reorder the second matrix to match the first one.
#' - \texttt{cossim} Vector of cosine similarities between the matched signatures.
#' - \texttt{distmat} Matrix of cosine similarities between all signatures.
#' 
#' @export
cosMatch <- function(H1,H2){
  if (!all.equal(dim(H1),dim(H2))){
    stop("The two signature matrices need to have the same dimensions")
  }
  K <- nrow(H1)
  d <- numeric(K)
  m <- numeric(K)
  dist <- sapply(1:K, function(y) sapply(1:K,function(x) similarity(H1[x,],H2[y,])))
  dist <- as.matrix(dist)
  distmat = dist
  for(s in 1:K){
    max.dist <- max(dist)
    remove = which(dist == max.dist, arr.ind = TRUE)
    dist[remove[1,1],] <- 0
    dist[,remove[1,2]] <- 0
    d[remove[1,1]] <- max.dist
    m[remove[1,1]] <- remove[1,2]
  }
  
  Output <- list()
  Output$match <- m          # the best matched signatures
  Output$cossim <- d          # Individual cosine similarity numbered after signatures in H1
  Output$distmat <- distmat 
  return(Output)
}