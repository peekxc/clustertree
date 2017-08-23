#include <Rcpp.h>
using namespace Rcpp;

#include <ANN/ANN.h> // ANN definitions
#include <ANN/kd_tree/kd_pr_search.h> // ANN kd tree priority search

List kdtree(NumericMatrix x, const int bkt_size = 30);
List kd_knn(NumericMatrix query_x, SEXP tree_ptr, int k, bool priority = false);