#include <Rcpp.h>
using namespace Rcpp;

List kdtree(NumericMatrix x);
List kd_knn(NumericMatrix query_x, SEXP tree_ptr, int k);