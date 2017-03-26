#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include <math.h>
#include "LSE.h"

using namespace Rcpp;

/**** algorithme forward pour le calcul de la vraisemblance ****/

// [[Rcpp::export]]
double logLikelihood(NumericMatrix logEmiss, NumericVector Dist, double a, double f) {
  double logf   = log(f);
  double logumf = log1p(-f);
  double lt00, lt01, lt10, lt11; // log proba transition
  int N = logEmiss.ncol();
  // état initial
  double alpha0 = logumf;
  double alpha1 = logf;
  double alpha0_;
  // Rcout << "alpha0 = " << alpha0 << ", alpha1 = " << alpha1 << "\n";
  if(f == 0) {
    for(int n = 0; n < N; n++) alpha0 += logEmiss(0,n);
    return alpha0;
  }
  if(f == 1) {
    for(int n = 0; n < N; n++) alpha1 += logEmiss(1,n);
    return alpha1;
  } 
  for(int n = 1; n < N; n++) {
    // calcul des lt
    double d = Dist[n-1];
    if(d < 0) { // changement de chromosomoe
      lt00 = logumf;
      lt01 = logf;
      lt10 = logumf;
      lt11 = logf;
    } else {
      // attention à la façon de calculer log(1 - exp(-a*d))
      // quand a petit, log1p(-exp(-a*d)) marche moins bien que log(-expm1(-a*d))
      // [enfin ça change pas grand chose]
      double ex = expm1(-a*d);
      lt00 = log1p(f*ex);
      lt01 = logf + log(-ex);
      lt10 = logumf + log(-ex);
      lt11 = log1p((1-f)*ex);
    }
    // calcul des alpha...
    alpha0_ = LSE( lt00 + logEmiss(0,n-1) + alpha0, lt10 + logEmiss(1,n-1) + alpha1 );
    alpha1  = LSE( lt01 + logEmiss(0,n-1) + alpha0, lt11 + logEmiss(1,n-1) + alpha1 );
    alpha0  = alpha0_;
  }
  return LSE(alpha0 + logEmiss(0,N-1), alpha1 + logEmiss(1,N-1));
}

RcppExport SEXP festim_logLikelihood(SEXP logEmissSEXP, SEXP DistSEXP, SEXP aSEXP, SEXP fSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type logEmiss(logEmissSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Dist(DistSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type f(fSEXP);
    __result = Rcpp::wrap(logLikelihood(logEmiss, Dist, a, f));
    return __result;
END_RCPP
}

