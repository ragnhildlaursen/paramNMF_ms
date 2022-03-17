#include <RcppArmadillo.h>
#include <cmath>        // std::abs
#include <tuple>
#include <iostream>


using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

double gkldev(arma::colvec y, arma::colvec mu) {
  int ySize = y.size();
  double sum = 0;
  for (int i=0; i<ySize; i++) {
    if (y[i] <= 0 || mu[i] <= 0) {
      sum += mu[i];
    } else {
      sum += y[i] * (log(y[i]) - log(mu[i])) - y[i] + mu[i];
    }
  }
  return sum;
}

// [[Rcpp::export]]
arma::rowvec glmUpdate(arma::colvec y, arma::mat X, int maxIter = 50, double epsilon=1e-8) {
  if (y.size() <= X.n_cols) {
    //Rcout << "Iterations:";
    return y.as_row();
  }
  //arma::colvec noise = arma::randu(y.size());
  arma::colvec mu = y + 0.1;
  double old = gkldev(y, mu);
  for (int i=0; i<maxIter; i++) {
    arma::colvec z = arma::log(mu) + (y-mu) / mu;
    arma::colvec w = arma::sqrt(mu);
    arma::colvec coef = arma::solve(X.each_col() % w, z % w);
    arma::colvec eta = X * coef;
    mu = arma::exp(eta);
    double newDev = gkldev(y, mu);
    if (2 * std::abs(old - newDev)/(0.1 + std::abs(2*newDev)) < epsilon) {
      break;
    }
    old = newDev;
  }
  return mu.as_row();
}

std::tuple<arma::mat, arma::mat, double> nmf1glm(arma::mat data, arma::mat designMatrices[], int noSignatures, int iter = 5000) {
  int genomes = data.n_rows;
  int mutTypes = data.n_cols;
  
  arma::mat exposures(genomes, noSignatures, arma::fill::randu);
  arma::mat signatures(noSignatures, mutTypes, arma::fill::randu);
  
  arma::mat estimate = exposures * signatures;
  
  for(int t = 0; t < iter; t++){
    exposures = exposures % ((data/estimate) * arma::trans(signatures));
    exposures = arma::normalise(exposures, 1);
    estimate = exposures * signatures;

    signatures = signatures % (arma::trans(exposures) * (data/estimate));

    for(int row = 0; row<noSignatures; row++) {
      signatures.row(row) = glmUpdate(arma::trans(signatures.row(row)), designMatrices[row]);
    }

    estimate = exposures * signatures;
  }
  double gkl = gkldev(arma::vectorise(data),arma::vectorise(estimate));
  
  return {exposures, signatures, gkl};
}

// [[Rcpp::export]]
List nmfprm(arma::mat data, List designMatrices, int noSignatures, int maxiter = 10000, double tolerance = 1e-8, int initial = 100, int smallIter = 500) {
  
  arma::mat designMatricesArray[noSignatures];
  for(int i=0; i<noSignatures; i++) {
    NumericMatrix x = designMatrices[i];
    designMatricesArray[i] = arma::mat(x.begin(), x.nrow(), x.ncol(), false);
  }  
  
  auto res = nmf1glm(data, designMatricesArray, noSignatures, smallIter);
  auto exposures = std::get<0>(res);
  auto signatures = std::get<1>(res);
  auto gklValue = std::get<2>(res);
  
  for(int i = 1; i < initial; i++){
    auto res = nmf1glm(data, designMatricesArray, noSignatures, smallIter);
    auto gklNew = std::get<2>(res);
    
    if(gklNew < gklValue){
      gklValue = gklNew;
      exposures = std::get<0>(res);
      signatures = std::get<1>(res);
    }
  }
  
  arma::mat estimate = exposures * signatures;
  
  double gklOld = gkldev(arma::vectorise(data),arma::vectorise(estimate));
  double gklNew = 2*gklOld;
  
  for(int t = 0; t < maxiter; t++){
    exposures = exposures % ((data/estimate) * arma::trans(signatures));
    exposures = arma::normalise(exposures, 1);
    estimate = exposures * signatures;
    
    
    signatures = signatures % (arma::trans(exposures) * (data/estimate));
    
    for(int row = 0; row<noSignatures; row++) {
      signatures.row(row) = glmUpdate(arma::trans(signatures.row(row)), designMatrices[row]);
    }
    
    estimate = exposures * signatures;
    
    gklNew = gkldev(arma::vectorise(data),arma::vectorise(estimate));
    
    if (2*std::abs(gklOld - gklNew)/(0.1 + std::abs(2*gklNew)) < tolerance){
      Rcout << "Iterations:";
      Rcout << t;
      Rcout << "\n";
      break;
    }
    gklOld = gklNew;
  }
  
  arma::colvec rsum = sum(signatures,1);
  exposures = exposures.each_row() % arma::trans(rsum);
  signatures = signatures.each_col() / rsum;
  
  
  List output = List::create(Named("exposures") = exposures,
                             Named("signatures") = signatures,
                             Named("gkl") = gklNew);
  return output;
}
