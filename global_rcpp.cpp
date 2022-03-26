#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <cmath>
#include <algorithm>
#include <vector>    
using namespace Rcpp;
using namespace arma;
// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//calculate Rt
// [[Rcpp::export]]
vec Rt_c(vec Rt, vec b, mat X_sd, int end, int n, int delta){
  for(int t=delta; t<end; t++){
    for(int i=0;i<n; i++){
      Rt[t]+=b[i]*X_sd(t-delta,i);
    }
    Rt[t]=exp(Rt[t]);
  }
  return Rt;
}


//calcualte logposterior
// [[Rcpp::export]]
double loglikelihood_c(int start, int end,double phi, double tau, vec Rt, vec w, vec I, vec I_imp){
  double loglikelihood=0;
  int n=static_cast<int>(w.size());
  for(int i=start; i<=end; i++){
    double m=0;
    for(int t=std::max(i-n,0);t<=i-1; t++){
      m+=I[t-1]*w[i-t-1]*Rt[t-1]+I_imp[t-1]*w[i-t-1]*phi;
    }
    if (m<0){
      return -99999;
    }
    loglikelihood+=R::dnbinom(I[i-1],m/(tau-1),1/tau,1);
  }
  return loglikelihood;
}