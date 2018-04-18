#ifndef ANN_UTIL_H
#define ANN_UTIL_H

#include "ANN.h"
#include "RcppHeader.h"

NumericVector ptToVec(ANNpoint pt, int dim);
IntegerVector idxArrayToVec(ANNidxArray idx_array, const int npts);
ANNpointArray matrixToANNpointArray(const NumericMatrix& x);

#endif