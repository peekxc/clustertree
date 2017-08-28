#ifndef R_KNN_H
#define R_KNN_H

#include "Rcpp.h"
#include "ANNx.h"

using namespace Rcpp;

static NumericMatrix emptyMatrix = NumericMatrix(0, 0);

// returns knn + dist
// [[Rcpp::export]]
List kNN_int(const NumericMatrix& x, int k,
             int type, int bucketSize, int splitRule, double approx,
             NumericMatrix& r_x); // put default in impl.

#endif