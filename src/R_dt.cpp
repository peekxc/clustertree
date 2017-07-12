#include <Rcpp.h>
using namespace Rcpp;

#include "dt_knn.h"

ANNpointArray matrixToANNpointArray(NumericMatrix x){
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

// R-facing API for the KNN-focused dual tree traversal
// [[Rcpp::export]]
List dt_knn(NumericMatrix q_x, const int k, NumericMatrix r_x = NumericMatrix(), const int bkt_size = 30, bool prune = false) {

  // If only query points are given, or q_x and r_x point to the same memory,
  // only one tree needs to be constructed
  bool identical_qr = r_x.size() <= 1 ? true : (&q_x) == (&r_x);
  // Rcout << "identical query and reference sets? " << identical_qr << " Size: " << r_x.size() << std::endl;
  ANNkd_tree* kd_treeQ, *kd_treeR;

  // Copy data over to ANN point array
  ANNpointArray qx_ann = matrixToANNpointArray(q_x);

  // Construct the dual tree KNN instance
  DualTreeKNN dt_knn = DualTreeKNN(prune);

  // Construct the tree(s)
  if (identical_qr){
    r_x = q_x; // Ensure r_x and q_x are identical incase r_x was NULL
    kd_treeQ = dt_knn.ConstructTree(qx_ann, q_x.nrow(), q_x.ncol(), bkt_size, ANN_KD_SUGGEST);
    kd_treeR = kd_treeQ;
  } else {
    ANNpointArray rx_ann = matrixToANNpointArray(r_x);
    kd_treeQ = dt_knn.ConstructTree(qx_ann, q_x.nrow(), q_x.ncol(), bkt_size, ANN_KD_SUGGEST);
    kd_treeR = dt_knn.ConstructTree(rx_ann, r_x.nrow(), r_x.ncol(), bkt_size, ANN_KD_SUGGEST);
  }

  // With the tree(s) created, setup KNN-specific bounds
  dt_knn.setup(kd_treeQ, kd_treeR);

  // Note: the search also returns the point itself (as the first hit)!
  // So we have to look for k+1 points.
  NumericMatrix dists(q_x.nrow(), k + 1);   // Distance matrix of kNN distances
  IntegerMatrix id(q_x.nrow(), k + 1);  // Id matrix of knn indices

  // Create dual tree using both trees
  dt_knn.KNN(k + 1, dists, id);
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

## Test for correctness
clustertree:::dt_knn(test_set, k = 4L, bkt_size = 1L, prune = FALSE)$dist[, -1] == dbscan::kNNdist(test_set, k = 4L)
clustertree:::dt_knn(test_set, k = 4L, bkt_size = 1L, prune = TRUE)$dist[, -1] == dbscan::kNNdist(test_set, k = 4L)
dbscan::kNN(test_set, k = 4L, sort = F, bucketSize = 1L)
## Benchmarking
size <- 1500
xyz <- as.matrix(data.frame(x = rnorm(size), y = rnorm(size), z = rnorm(size)))
microbenchmark::microbenchmark(invisible(clustertree:::dt_knn(xyz, k = 30L, bkt_size = 15L, prune = FALSE)), times = 15L)
microbenchmark::microbenchmark(invisible(clustertree:::dt_knn(xyz, k = 30L, bkt_size = 15L, prune = TRUE)), times = 15L)
microbenchmark::microbenchmark(invisible(dbscan::kNN(xyz, k = 30L, sort = F, bucketSize = 15L)), times = 15L)

*/
