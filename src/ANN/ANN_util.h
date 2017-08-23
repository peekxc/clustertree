#ifndef ANN_UTIL_H
#define ANN_UTIL_H

#include <ANN/ANN.h>
#include <Rcpp.h>
using namespace Rcpp;

static ANNpointArray matrixToANNpointArray(NumericMatrix x){
  // Copy data over to ANN point array
  int nrow = x.nrow(), ncol = x.ncol();
  ANNpointArray dataPts = annAllocPts(nrow, ncol);
  for(int i = 0; i < nrow; i++){
    for(int j = 0; j < ncol; j++){
      (dataPts[i])[j] = x(i, j);
    }
  }
  return dataPts;
}

#endif