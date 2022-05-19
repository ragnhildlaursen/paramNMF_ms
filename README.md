# paramNMF_ms
This project parametrizes mutational signatures to give more robust and interpretable results. 

The data and codes here can be used to recreate the results seen in the paper 'Robust estimation of mutational signatures by adaptive learning of nucleotide interaction terms in non–negative matrix factorization' by Ragnhild Laursen, Lasse Marretty and Asger Hobolth. 

The code files **ModelSelection.R** and **NMF2.cpp** are the essential code used to run non-negative matrix factorization with parametrization of the signatures. The ModelSelection.R file include the original codes and NMF2.cpp includes a main function **nmfprm()** to find the estimates, where it is implemented i C++ to achieve faster convergence.  

## Analysis of BRCA datasets 
The code for recreating the results for the BRCA data sets can be found in the *BRCA* folder. This folder consist of both data sets, codes to create the results and a folder that holds the results. 

The two breast cancer datasets of mutational counts are **BRCA21.RData** and **BRCA214.RData**. They consist of a datamatrix of size $21 \times 96$ and $214 \times 96$, respectively.

The coding files do the folowing:
 - **loadBRCA21models.R** - loads the BRCA21 dataset and creates the factors and model matrices for the parametrization of the signatures.
 - **loadBRCA214models.R** - loads the BRCA214 dataset and creates the factors and model matrices for the parametrization of the signatures.
 - **optimalSignaturesBRCA.R** - recovers the fit and estimates of $W$ and $H$ for each of the parametrized models. 
 - **RobustExposures.R** - find the cosine similarity between the orignal exposures and the reestimated ones from the downsampled data. 
 - **RobustSignatures.R** - find the cosine similarity between the original signatures and the reestimed ones after applying parametric bootstrapping to the original results.

## Analysis of UCUT datasets

The folder *UCUT* contains all the necessary data and code to recreate the results for the UCUT analysis. The structure of the folder is equivalent to the *BRCA* folder, where the dataset is the file **UCUT_5_all.RData**. The data matrix here is of size $26 \times 1536$.


## Code for plotting the figures 
The code for plotting the figures in the paper are found in *plotcode* and the resulting plots can be found in the folder *plots*.
