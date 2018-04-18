#include "ANN_util.h"

NumericVector ptToVec(ANNpoint pt, int dim){
  NumericVector new_pt = NumericVector(dim);
  for (int d_i = 0; d_i < dim; ++d_i){
    new_pt[d_i] = pt[d_i];
  }
  return(new_pt);
}

IntegerVector idxArrayToVec(ANNidxArray idx_array, const int npts){
  IntegerVector idx_vec = IntegerVector(npts);
  for (int i = 0; i < npts; ++i){ idx_vec[i] = idx_array[i]; }
  return(idx_vec);
}

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