// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Rginv
arma::mat Rginv(const arma::mat& m);
RcppExport SEXP _kalmanfilter_Rginv(SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(Rginv(m));
    return rcpp_result_gen;
END_RCPP
}
// gen_inv
arma::mat gen_inv(arma::mat& m);
RcppExport SEXP _kalmanfilter_gen_inv(SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(gen_inv(m));
    return rcpp_result_gen;
END_RCPP
}
// contains
bool contains(std::string s, Rcpp::List L);
RcppExport SEXP _kalmanfilter_contains(SEXP sSEXP, SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type s(sSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type L(LSEXP);
    rcpp_result_gen = Rcpp::wrap(contains(s, L));
    return rcpp_result_gen;
END_RCPP
}
// kalman_filter_cpp
Rcpp::List kalman_filter_cpp(Rcpp::List& ssm, const arma::mat& yt, Rcpp::Nullable<Rcpp::NumericMatrix> Xo, Rcpp::Nullable<Rcpp::NumericMatrix> Xs, Rcpp::Nullable<Rcpp::NumericMatrix> weight, bool smooth);
RcppExport SEXP _kalmanfilter_kalman_filter_cpp(SEXP ssmSEXP, SEXP ytSEXP, SEXP XoSEXP, SEXP XsSEXP, SEXP weightSEXP, SEXP smoothSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List& >::type ssm(ssmSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type yt(ytSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericMatrix> >::type Xo(XoSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericMatrix> >::type Xs(XsSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericMatrix> >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< bool >::type smooth(smoothSEXP);
    rcpp_result_gen = Rcpp::wrap(kalman_filter_cpp(ssm, yt, Xo, Xs, weight, smooth));
    return rcpp_result_gen;
END_RCPP
}
// kalman_filter_tvp_cpp
Rcpp::List kalman_filter_tvp_cpp(Rcpp::List& ssm, const arma::mat& yt, Rcpp::Nullable<Rcpp::NumericMatrix> Xo, Rcpp::Nullable<Rcpp::NumericMatrix> Xs, Rcpp::Nullable<Rcpp::NumericMatrix> weight, bool smooth);
RcppExport SEXP _kalmanfilter_kalman_filter_tvp_cpp(SEXP ssmSEXP, SEXP ytSEXP, SEXP XoSEXP, SEXP XsSEXP, SEXP weightSEXP, SEXP smoothSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List& >::type ssm(ssmSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type yt(ytSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericMatrix> >::type Xo(XoSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericMatrix> >::type Xs(XsSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericMatrix> >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< bool >::type smooth(smoothSEXP);
    rcpp_result_gen = Rcpp::wrap(kalman_filter_tvp_cpp(ssm, yt, Xo, Xs, weight, smooth));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_kalmanfilter_Rginv", (DL_FUNC) &_kalmanfilter_Rginv, 1},
    {"_kalmanfilter_gen_inv", (DL_FUNC) &_kalmanfilter_gen_inv, 1},
    {"_kalmanfilter_contains", (DL_FUNC) &_kalmanfilter_contains, 2},
    {"_kalmanfilter_kalman_filter_cpp", (DL_FUNC) &_kalmanfilter_kalman_filter_cpp, 6},
    {"_kalmanfilter_kalman_filter_tvp_cpp", (DL_FUNC) &_kalmanfilter_kalman_filter_tvp_cpp, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_kalmanfilter(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
