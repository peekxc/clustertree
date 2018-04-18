#include "RcppHeader.h"
#include "ANN.h" // ANN definitions
#include "kd_pr_search.h" // ANN kd tree priority search

List kdtree(NumericMatrix x, const int bkt_size = 30);
List kd_knn(NumericMatrix query_x, SEXP tree_ptr, int k, bool priority = false);