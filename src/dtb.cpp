#include <Rcpp.h>
using namespace Rcpp;

#include <iostream>     // std::cout, std::ostream, std::ios
#include <fstream>      // std::filebuf

// ANN library
#include "ANN/ANN.h"
#include "kd_search.h" // kd-search declarations
#include "bd_tree.h"	 // bd-tree declarations

// KD tree R-accessible functions
#include "R_kdtree.h"

// ANN extensions
#include "dual_tree.h"

// [[Rcpp::export]]
List DT_knn(NumericMatrix x, const int k) {

  // Copy data over to ANN point array
  int nrow = x.nrow(), ncol = x.ncol();

  ANNpointArray dataPts = annAllocPts(nrow, ncol);
  for(int i = 0; i < nrow; i++){
    for(int j = 0; j < ncol; j++){
      (dataPts[i])[j] = x(i, j);
    }
  }

  // Create kd tree
  ANNkd_tree* kdTree = new ANNkd_tree(dataPts, nrow, ncol, 30, ANN_KD_SUGGEST);

  // Note: the search also returns the point itself (as the first hit)!
  // So we have to look for k+1 points.
  // ANNdistArray dists = new ANNdist[k+1];
  // ANNidxArray nnIdx = new ANNidx[k+1];

  // Distance matrix of kNN distances
  NumericMatrix dists(x.nrow(), k);

  // Id matrix of knn indices
  IntegerMatrix id(x.nrow(), k);

  // Create dual tree using both trees
  DualTree dt = DualTree(kdTree, kdTree, k);
  List test_status = dt.test_cases();
  // dt.KNN(k, dists, id);
  List res = List::create(_["dist"] = dists, _["id"] = id, _["info"] = test_status);

  return(res);
}

/*** R
  test_set <- matrix(c(13, 291,
                       57, 145,
                       115, 232,
                       86, 27,
                       145, 28,
                       262, 203,
                       320, 261,
                       379, 174,
                       261, 71,
                       325, 57), byrow=T, ncol=2)
  test_set[, 2] <- 321 - test_set[, 2]
  clustertree:::DT_knn(test_set, k = 3L)
*/
