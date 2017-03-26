#include <Rcpp.h>
using namespace Rcpp;

// i et j sont des états 0 / 1
// proba d'une transition j -> i
double transition_proba(int j, int i, double d, double a, double f) {
  if(d < 0) {
    if(i)
      return(f);
    return(1-f);
  }
  if(i) {
    if(j)
      return(1 + (1-f)*expm1(-a*d)); // t11
    return(-f*expm1(-a*d));    // t01
  }
  if(j)
    return(-(1-f)*expm1(-a*d)); // t10
  return(1 + f*expm1(-a*d));    // t00
}

// On ne simule que la chaîne S
NumericVector simu0(NumericVector Dist, double a, double f) {
  int N = Dist.length()+1;
  NumericVector S(N);
  int state = S[0] = (R::runif(0,1)<f)?1:0;
  for(int n = 1; n < N; n++) {
    double tr = transition_proba(state,0,Dist[n-1],a,f);
    state = S[n] = R::runif(0,1)<tr?0:1;
  }
  return S;
}

// ** étant donné l'état, simuler les génotypes... **

// ici, juste un allèle
int simu1a(NumericMatrix Freq, int index) {
  double s = 0;
  int M = Freq.ncol();
  double r = R::runif(0,1);
  for(int i = 0; i < M; i++) {
    s += Freq(index, i);
    if(s >= r) return i+1;
  }
  return M; // Ou message d'erreur... 
}

// [[Rcpp::export]]
NumericMatrix simu_geno(NumericVector S, NumericMatrix Freq) {
  int N = S.length();
  if(N != Freq.nrow()) stop("Dimensions mismatch");
  NumericMatrix Y(N,2);
  for(int i = 0; i < N; i++) {
    Y(i,0) = simu1a(Freq, i);
    if(S(i) == 0) 
      Y(i,1) =  simu1a(Freq, i);
    else 
      Y(i,1) = Y(i,0);
  }
  return Y;
}


RcppExport SEXP festim_simu0(SEXP DistSEXP, SEXP aSEXP, SEXP fSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type Dist(DistSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type f(fSEXP);
    __result = Rcpp::wrap(simu0(Dist, a, f));
    return __result;
END_RCPP
}

RcppExport SEXP festim_simu_geno(SEXP SSEXP, SEXP FreqSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type S(SSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Freq(FreqSEXP);
    __result = Rcpp::wrap(simu_geno(S, Freq));
    return __result;
END_RCPP
}


