// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;


RcppExport SEXP _rcpp_module_boot_stan_fit4linear_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4splines_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4splines_interaction_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4splines_interaction_sparse_mod();

static const R_CallMethodDef CallEntries[] = {
    {"_rcpp_module_boot_stan_fit4linear_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4linear_mod, 0},
    {"_rcpp_module_boot_stan_fit4splines_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4splines_mod, 0},
    {"_rcpp_module_boot_stan_fit4splines_interaction_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4splines_interaction_mod, 0},
    {"_rcpp_module_boot_stan_fit4splines_interaction_sparse_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4splines_interaction_sparse_mod, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_treeRingSplines(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
