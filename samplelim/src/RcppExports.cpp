// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// inner_ball
Rcpp::NumericVector inner_ball(Rcpp::Reference P);
RcppExport SEXP _samplelim_inner_ball(SEXP PSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Reference >::type P(PSEXP);
    rcpp_result_gen = Rcpp::wrap(inner_ball(P));
    return rcpp_result_gen;
END_RCPP
}
// sample_points
Rcpp::NumericMatrix sample_points(Rcpp::Reference P, int n, Rcpp::Nullable<Rcpp::List> random_walk);
RcppExport SEXP _samplelim_sample_points(SEXP PSEXP, SEXP nSEXP, SEXP random_walkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Reference >::type P(PSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type random_walk(random_walkSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_points(P, n, random_walk));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_samplelim_inner_ball", (DL_FUNC) &_samplelim_inner_ball, 1},
    {"_samplelim_sample_points", (DL_FUNC) &_samplelim_sample_points, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_samplelim(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
