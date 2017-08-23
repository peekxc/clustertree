// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// kruskalsMST
NumericMatrix kruskalsMST(const NumericVector dist_x);
RcppExport SEXP _clustertree_kruskalsMST(SEXP dist_xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type dist_x(dist_xSEXP);
    rcpp_result_gen = Rcpp::wrap(kruskalsMST(dist_x));
    return rcpp_result_gen;
END_RCPP
}
// primsMST
NumericMatrix primsMST(const NumericVector dist_x);
RcppExport SEXP _clustertree_primsMST(SEXP dist_xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type dist_x(dist_xSEXP);
    rcpp_result_gen = Rcpp::wrap(primsMST(dist_x));
    return rcpp_result_gen;
END_RCPP
}
// primsRSL
NumericMatrix primsRSL(const NumericVector r, const NumericVector r_k, const int n, const double alpha, const int type);
RcppExport SEXP _clustertree_primsRSL(SEXP rSEXP, SEXP r_kSEXP, SEXP nSEXP, SEXP alphaSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type r_k(r_kSEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const int >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(primsRSL(r, r_k, n, alpha, type));
    return rcpp_result_gen;
END_RCPP
}
// mstToHclust
List mstToHclust(NumericMatrix mst, const int n);
RcppExport SEXP _clustertree_mstToHclust(SEXP mstSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type mst(mstSEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(mstToHclust(mst, n));
    return rcpp_result_gen;
END_RCPP
}
// cut_simplified_hclust
IntegerVector cut_simplified_hclust(List hcl, IntegerVector cl_in, const int big_n);
RcppExport SEXP _clustertree_cut_simplified_hclust(SEXP hclSEXP, SEXP cl_inSEXP, SEXP big_nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type hcl(hclSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type cl_in(cl_inSEXP);
    Rcpp::traits::input_parameter< const int >::type big_n(big_nSEXP);
    rcpp_result_gen = Rcpp::wrap(cut_simplified_hclust(hcl, cl_in, big_n));
    return rcpp_result_gen;
END_RCPP
}
// simplified_hclust
List simplified_hclust(List hcl, const int min_sz);
RcppExport SEXP _clustertree_simplified_hclust(SEXP hclSEXP, SEXP min_szSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type hcl(hclSEXP);
    Rcpp::traits::input_parameter< const int >::type min_sz(min_szSEXP);
    rcpp_result_gen = Rcpp::wrap(simplified_hclust(hcl, min_sz));
    return rcpp_result_gen;
END_RCPP
}
// clusterTree
List clusterTree(const NumericVector dist_x, const NumericVector r_k, const int k, const double alpha, const int type);
RcppExport SEXP _clustertree_clusterTree(SEXP dist_xSEXP, SEXP r_kSEXP, SEXP kSEXP, SEXP alphaSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type dist_x(dist_xSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type r_k(r_kSEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const int >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(clusterTree(dist_x, r_k, k, alpha, type));
    return rcpp_result_gen;
END_RCPP
}
// naive_clustertree
List naive_clustertree(const NumericVector x, const NumericVector r_k, const double alpha, const int type);
RcppExport SEXP _clustertree_naive_clustertree(SEXP xSEXP, SEXP r_kSEXP, SEXP alphaSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type r_k(r_kSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const int >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(naive_clustertree(x, r_k, alpha, type));
    return rcpp_result_gen;
END_RCPP
}
// pruneCT
NumericVector pruneCT(List C_n, NumericVector prune_heights, IntegerVector valid_idx);
RcppExport SEXP _clustertree_pruneCT(SEXP C_nSEXP, SEXP prune_heightsSEXP, SEXP valid_idxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type C_n(C_nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prune_heights(prune_heightsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type valid_idx(valid_idxSEXP);
    rcpp_result_gen = Rcpp::wrap(pruneCT(C_n, prune_heights, valid_idx));
    return rcpp_result_gen;
END_RCPP
}
// vol_nSphere
double vol_nSphere(const int n, const double R);
RcppExport SEXP _clustertree_vol_nSphere(SEXP nSEXP, SEXP RSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const double >::type R(RSEXP);
    rcpp_result_gen = Rcpp::wrap(vol_nSphere(n, R));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_clustertree_kruskalsMST", (DL_FUNC) &_clustertree_kruskalsMST, 1},
    {"_clustertree_primsMST", (DL_FUNC) &_clustertree_primsMST, 1},
    {"_clustertree_primsRSL", (DL_FUNC) &_clustertree_primsRSL, 5},
    {"_clustertree_mstToHclust", (DL_FUNC) &_clustertree_mstToHclust, 2},
    {"_clustertree_cut_simplified_hclust", (DL_FUNC) &_clustertree_cut_simplified_hclust, 3},
    {"_clustertree_simplified_hclust", (DL_FUNC) &_clustertree_simplified_hclust, 2},
    {"_clustertree_clusterTree", (DL_FUNC) &_clustertree_clusterTree, 5},
    {"_clustertree_naive_clustertree", (DL_FUNC) &_clustertree_naive_clustertree, 4},
    {"_clustertree_pruneCT", (DL_FUNC) &_clustertree_pruneCT, 3},
    {"_clustertree_vol_nSphere", (DL_FUNC) &_clustertree_vol_nSphere, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_clustertree(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
