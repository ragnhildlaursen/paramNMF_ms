
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
double gkldev(arma::vec y, arma::vec mu){
  double sum = 0;
  int ny = y.size();
  for(int i = 0; i < ny; i++){
    if (y[i] > 0){
      sum +=  y[i] * log(y[i]/mu[i]) - y[i] + mu[i];
    }else{
      sum += mu[i];
    }
  }
  return sum;
}



// [[Rcpp::export]]
List nmf1(arma::mat data, int noSignatures, int iter = 5000) {
  int genomes = data.n_rows;
  int mutTypes = data.n_cols;
  
  arma::mat Exposures(genomes, noSignatures, arma::fill::randu);
  arma::mat Signatures(noSignatures, mutTypes, arma::fill::randu);
  
  arma::mat estimate = Exposures * Signatures;
  
  for(int t = 0; t < iter; t ++){
    Exposures = Exposures % ((data/estimate) * trans(Signatures));
    estimate = Exposures * Signatures;
    
    Signatures = Signatures % (trans(Exposures) * (data/estimate));
    Signatures = normalise(Signatures, 1, 1);
    estimate = Exposures * Signatures;
  }
  double gkl = gkldev(arma::vectorise(data),arma::vectorise(estimate));
  
  List Output = List::create(Named("Exposures") = Exposures,
                             Named("Signatures") = Signatures,
                             Named("gkl") = gkl);
  return Output;
}


// [[Rcpp::export]]
List nmf2(arma::mat data, int noSignatures, int maxiter = 5000, double tolerance = 1e-8, int initial = 100) {
  
  List res = nmf1(data, noSignatures, 100);
  
  double gklvalue = res["gkl"];
  
  arma::mat Exposures = res["Exposures"];
  arma::mat Signatures = res["Signatures"];
  
  for(int i = 1; i < initial; i++){
    res = nmf1(data, noSignatures, 100);
    double gklnew = res["gkl"];
    
    if(gklnew < gklvalue){
      gklvalue = res["gkl"];
      arma::mat Exposures = res["Exposures"];
      arma::mat Signatures = res["Signatures"];
    }
  }
  
  arma::mat estimate = Exposures * Signatures;
  
  double gklold = gkldev(arma::vectorise(data),arma::vectorise(estimate));
  double gklnew = 2*gklold;
  
  for(int t = 0; t < maxiter; t ++){
    Exposures = Exposures % ((data/estimate) * trans(Signatures));
    estimate = Exposures * Signatures;
    
    Signatures = Signatures % (trans(Exposures) * (data/estimate));
    Signatures = normalise(Signatures, 1, 1);
    estimate = Exposures * Signatures;
    
    gklnew = gkldev(arma::vectorise(data),arma::vectorise(estimate));
    
    if(2*std::abs(gklold - gklnew)/(0.1 + std::abs(2*gklnew)) < tolerance) break;
    gklold = gklnew;
  }
  
  List Output = List::create(Named("Exposures") = Exposures,
                             Named("Signatures") = Signatures,
                             Named("gkl") = gklnew);
  return Output;
}
// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
res = nmf1(matrix(runif(9), nrow = 3), noSig = 2)
rowSums(res$Signatures)
nmf2(matrix(runif(9), nrow = 3), noSig = 2)
*/
