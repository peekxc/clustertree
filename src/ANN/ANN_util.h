#ifndef ANN_UTIL_H
#define ANN_UTIL_H

#include <ANN/ANN.h>
#include <Rcpp.h>
using namespace Rcpp;

ANNpointArray matrixToANNpointArray(const NumericMatrix& x);

#endif