#include <Rcpp.h>
using namespace Rcpp;

#include "dt_knn.h"

// R-facing API for the KNN-focused dual tree traversal
// [[Rcpp::export]]
List dt_knn(NumericMatrix x, const int k, const int bkt_size = 30, bool prune = false) {

  // Copy data over to ANN point array
  int nrow = x.nrow(), ncol = x.ncol();
  ANNpointArray dataPts = annAllocPts(nrow, ncol);
  for(int i = 0; i < nrow; i++){
    for(int j = 0; j < ncol; j++){
      (dataPts[i])[j] = x(i, j);
    }
  }

  // Create kd tree
  ANNkd_tree* kdTree = new ANNkd_tree(dataPts, nrow, ncol, bkt_size, ANN_KD_SUGGEST);

  // Note: the search also returns the point itself (as the first hit)!
  // So we have to look for k+1 points.
  NumericMatrix dists(x.nrow(), k + 1);   // Distance matrix of kNN distances
  IntegerMatrix id(x.nrow(), k + 1);  // Id matrix of knn indices

  // Create dual tree using both trees
  DualTreeKNN dt_knn = DualTreeKNN(kdTree, kdTree);
  dt_knn.KNN(k + 1, dists, id, prune);
  List res = List::create(_["dist"] = dists, _["id"] = id + 1);


  // Performance test cases
  // // Test cases
  // List test_res = List();
  // dt.test_cases(test_res);
  // res["test_results] = test_res;


  // Return results
  return(res);
}

/*** R
## Testing the search
test_set <- matrix(c(13, 291, 57, 145, 115, 232, 86, 27, 145, 28, 262, 203, 320, 261, 379, 174, 261, 71, 325, 57), byrow=T, ncol=2)
test_set[, 2] <- 321 - test_set[, 2]
clustertree:::dt_knn(test_set, k = 4L, bkt_size = 1L, prune = FALSE)
clustertree:::dt_knn(test_set, k = 4L, bkt_size = 1L, prune = TRUE)


## Benchmarking
size <- 1500
xyz <- as.matrix(data.frame(x = rnorm(size), y = rnorm(size), z = rnorm(size)))
microbenchmark::microbenchmark(invisible(clustertree:::dt_knn(xyz, k = 30L, bkt_size = 15L, prune = FALSE)), times = 15L)
microbenchmark::microbenchmark(invisible(clustertree:::dt_knn(xyz, k = 30L, bkt_size = 15L, prune = TRUE)), times = 15L)
microbenchmark::microbenchmark(invisible(dbscan::kNN(xyz, k = 30L, sort = F, bucketSize = 15L)), times = 15L)

*/
