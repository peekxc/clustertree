#include "ANN_util.h"

ANNpointArray matrixToANNpointArray(const NumericMatrix& x){
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