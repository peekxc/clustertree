// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// kruskalsMST
NumericMatrix kruskalsMST(const NumericVector dist_x);
RcppExport SEXP clustertree_kruskalsMST(SEXP dist_xSEXP) {
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
RcppExport SEXP clustertree_primsMST(SEXP dist_xSEXP) {
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
RcppExport SEXP clustertree_primsRSL(SEXP rSEXP, SEXP r_kSEXP, SEXP nSEXP, SEXP alphaSEXP, SEXP typeSEXP) {
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
RcppExport SEXP clustertree_mstToHclust(SEXP mstSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type mst(mstSEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(mstToHclust(mst, n));
    return rcpp_result_gen;
END_RCPP
}
// clusterTree
List clusterTree(const NumericVector dist_x, const NumericVector r_k, const int k, const double alpha, const int type, IntegerVector knn_indices);
RcppExport SEXP clustertree_clusterTree(SEXP dist_xSEXP, SEXP r_kSEXP, SEXP kSEXP, SEXP alphaSEXP, SEXP typeSEXP, SEXP knn_indicesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type dist_x(dist_xSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type r_k(r_kSEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type knn_indices(knn_indicesSEXP);
    rcpp_result_gen = Rcpp::wrap(clusterTree(dist_x, r_k, k, alpha, type, knn_indices));
    return rcpp_result_gen;
END_RCPP
}
// DT_knn
List DT_knn(NumericMatrix x, const int k, const int bkt_size);
RcppExport SEXP clustertree_DT_knn(SEXP xSEXP, SEXP kSEXP, SEXP bkt_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    Rcpp::traits::input_parameter< const int >::type bkt_size(bkt_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(DT_knn(x, k, bkt_size));
    return rcpp_result_gen;
END_RCPP
}
// naive_clustertree
List naive_clustertree(const NumericVector x, const NumericVector r_k, const double alpha, const int type);
RcppExport SEXP clustertree_naive_clustertree(SEXP xSEXP, SEXP r_kSEXP, SEXP alphaSEXP, SEXP typeSEXP) {
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
// kd_knn
List kd_knn(NumericMatrix query_x, SEXP tree_ptr, int k);
RcppExport SEXP clustertree_kd_knn(SEXP query_xSEXP, SEXP tree_ptrSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type query_x(query_xSEXP);
    Rcpp::traits::input_parameter< SEXP >::type tree_ptr(tree_ptrSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(kd_knn(query_x, tree_ptr, k));
    return rcpp_result_gen;
END_RCPP
}
// kdtree
List kdtree(NumericMatrix x);
RcppExport SEXP clustertree_kdtree(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(kdtree(x));
    return rcpp_result_gen;
END_RCPP
}
// kNN_int
List kNN_int(NumericMatrix data, int k, int type, int bucketSize, int splitRule, double approx);
RcppExport SEXP clustertree_kNN_int(SEXP dataSEXP, SEXP kSEXP, SEXP typeSEXP, SEXP bucketSizeSEXP, SEXP splitRuleSEXP, SEXP approxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< int >::type bucketSize(bucketSizeSEXP);
    Rcpp::traits::input_parameter< int >::type splitRule(splitRuleSEXP);
    Rcpp::traits::input_parameter< double >::type approx(approxSEXP);
    rcpp_result_gen = Rcpp::wrap(kNN_int(data, k, type, bucketSize, splitRule, approx));
    return rcpp_result_gen;
END_RCPP
}
// test_ptr
IntegerVector test_ptr();
RcppExport SEXP clustertree_test_ptr() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(test_ptr());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"clustertree_kruskalsMST", (DL_FUNC) &clustertree_kruskalsMST, 1},
    {"clustertree_primsMST", (DL_FUNC) &clustertree_primsMST, 1},
    {"clustertree_primsRSL", (DL_FUNC) &clustertree_primsRSL, 5},
    {"clustertree_mstToHclust", (DL_FUNC) &clustertree_mstToHclust, 2},
    {"clustertree_clusterTree", (DL_FUNC) &clustertree_clusterTree, 6},
    {"clustertree_DT_knn", (DL_FUNC) &clustertree_DT_knn, 3},
    {"clustertree_naive_clustertree", (DL_FUNC) &clustertree_naive_clustertree, 4},
    {"clustertree_kd_knn", (DL_FUNC) &clustertree_kd_knn, 3},
    {"clustertree_kdtree", (DL_FUNC) &clustertree_kdtree, 1},
    {"clustertree_kNN_int", (DL_FUNC) &clustertree_kNN_int, 6},
    {"clustertree_test_ptr", (DL_FUNC) &clustertree_test_ptr, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_clustertree(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
