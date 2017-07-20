#include <Rcpp.h>
using namespace Rcpp;

List kdtree(NumericMatrix x, const int bkt_size = 30);
List kd_knn(NumericMatrix query_x, SEXP tree_ptr, int k, bool priority = false);