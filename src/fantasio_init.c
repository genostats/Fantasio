#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/*
  generated with tools::package_native_routine_registration_skeleton
*/

/* .Call calls */
extern SEXP festim_forward_backward(SEXP, SEXP, SEXP, SEXP);
extern SEXP festim_logEmiss(SEXP, SEXP, SEXP, SEXP);
extern SEXP festim_logLikelihood(SEXP, SEXP, SEXP, SEXP);
extern SEXP festim_logLikelihood_gradient(SEXP, SEXP, SEXP, SEXP);
extern SEXP festim_m4_logEmiss(SEXP, SEXP, SEXP, SEXP);
extern SEXP festim_simu0(SEXP, SEXP, SEXP);
extern SEXP festim_simu_geno(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"festim_forward_backward",       (DL_FUNC) &festim_forward_backward,       4},
    {"festim_logEmiss",               (DL_FUNC) &festim_logEmiss,               4},
    {"festim_logLikelihood",          (DL_FUNC) &festim_logLikelihood,          4},
    {"festim_logLikelihood_gradient", (DL_FUNC) &festim_logLikelihood_gradient, 4},
    {"festim_m4_logEmiss",            (DL_FUNC) &festim_m4_logEmiss,            4},
    {"festim_simu0",                  (DL_FUNC) &festim_simu0,                  3},
    {"festim_simu_geno",              (DL_FUNC) &festim_simu_geno,              2},
    {NULL, NULL, 0}
};

void R_init_Fantasio(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

