#include <Rcpp.h>
using namespace Rcpp;

#include <utilities.h> // R_INFO, profiling mode, etc. must include first
#include <DT/KNN/dt_knn.h> // Dual Tree KNN implementation
#include "R_kdtree.h" // R Interface to ANN KD trees
#include <ANN/ANN_util.h>

// R-facing API for the KNN-focused dual tree traversal
// [[Rcpp::export]]
List dt_knn(NumericMatrix q_x, const int k, NumericMatrix r_x = NumericMatrix(), const int bkt_size = 30, bool prune = true) {

  // If only query points are given, or q_x and r_x point to the same memory,
  // only one tree needs to be constructed
  bool identical_qr = r_x.size() <= 1 ? true : (&q_x) == (&r_x);
  ANNkd_tree* kd_treeQ, *kd_treeR;

  // Copy data over to ANN point array
  ANNpointArray qx_ann = matrixToANNpointArray(q_x);

  // Construct the dual tree KNN instance
  const int d = q_x.ncol();
  L_1 metric = L_1(d);
  DualTreeKNN dt_knn = DualTreeKNN(prune, d, metric);

  // Construct the tree(s)
  if (identical_qr){
    r_x = q_x; // Ensure r_x and q_x are identical incase r_x was NULL
    BEGIN_PROFILE()
    kd_treeQ = dt_knn.ConstructTree(qx_ann, q_x.nrow(), q_x.ncol(), bkt_size, ANN_KD_SUGGEST);
    kd_treeR = kd_treeQ;
    REPORT_TIME("KD Tree construction")
  } else {
    ANNpointArray rx_ann = matrixToANNpointArray(r_x);
    kd_treeQ = dt_knn.ConstructTree(qx_ann, q_x.nrow(), q_x.ncol(), bkt_size, ANN_KD_SUGGEST);
    kd_treeR = dt_knn.ConstructTree(rx_ann, r_x.nrow(), r_x.ncol(), bkt_size, ANN_KD_SUGGEST);
  }

  // With the tree(s) created, setup KNN-specific bounds
  dt_knn.setup(kd_treeQ, kd_treeR);

  // Note: the search also returns the point itself (as the first hit)!
  // So we have to look for k+1 points.

  // Create dual tree using both trees
  UTIL(dt_knn.PrintTree((ANNbool) true, true))

  // ----- Start KNN Performance Testing -----
  List res;
  #ifdef PROFILING
  annResetStats(r_x.nrow());

  // Regular kdtree search
  List ann_kdtree = kdtree(r_x, bkt_size); // create regular kdtree with the same bucket size
  BEGIN_PROFILE()
    List res1 = kd_knn(q_x, (SEXP) ann_kdtree["kdtree_ptr"], k, false);
  REPORT_TIME("Regular KNN search")

  // Print statistics
  Rcout << "Regular KDtree search performance: " << std::endl;
  annPrintStats((ANNbool) false);

  // Dual tree equivalent
  Rcout << "DualTree search performance: " << std::endl;
  annResetStats(r_x.nrow());
  annResetCounts();	// reset stats for a set of queries
  BEGIN_PROFILE()
  res = dt_knn.KNN(k+1);
  REPORT_TIME("Dual tree KNN search")
  annUpdateStats();
  annPrintStats((ANNbool) false);
  #endif
  // ----- End Performance Testing -----


  // Return results
  return(res);
}

/*** R
## Testing the search
test_set <- matrix(c(13, 291, 57, 145, 115, 232, 86, 27, 145, 28, 262, 203, 320, 261, 379, 174, 261, 71, 325, 57), byrow=T, ncol=2)
test_set[, 2] <- 321 - test_set[, 2]
clustertree:::dt_knn(test_set, k = 4L, bkt_size = 1L, prune = TRUE)
clustertree:::dt_knn(test_set, k = 4L, bkt_size = 1L, prune = FALSE)



test_knn <- clustertree:::dt_knn(ts2, k = 40L, bkt_size = 1L, prune = TRUE)


l1_dist <- as.matrix(dist(ts2, method = "manhattan"))
l1_knn_id <- t(apply(l1_dist, 1, order))
l1_knn_dist <- t(apply(l1_dist, 1, sort))

all(test_knn$dist == l1_knn_dist[, 2:41L])
all(test_knn$id == l1_knn_id[, 2:41L])

size <- 50
ts2 <- as.matrix(data.frame(x = rnorm(size), y = rnorm(size)))
invisible(clustertree:::dt_knn(ts2, k = 15L, bkt_size = 5L, prune = TRUE))

## Test for correctness
bkt_sz = 5L
clustertree:::dt_knn(test_set, k = 8L, bkt_size = bkt_sz, prune = TRUE)$dist == dbscan::kNNdist(test_set, k = 8L)
clustertree:::dt_knn(test_set, k = 8L, bkt_size = bkt_sz, prune = TRUE)$id == unname(dbscan::kNN(test_set, k = 8L)$id)


#sapply(seq(1, 140, by = 10), function(i){
  X_n <- as.matrix(iris[c(41, 42, 43, 44, 45, 46, 47, 50, 51), 1:4])
  k <- 2L
  truth <- dbscan::kNN(X_n, k = k)
  all(clustertree:::dt_knn(X_n, k = k, bkt_size = bkt_sz, prune = TRUE)$dist == truth$dist)
  #clustertree:::dt_knn(X_n, k = k, bkt_size = bkt_sz, prune = TRUE)$id == unname(truth$id)
#})


## Benchmarking
size <- 1500
xyz <- as.matrix(data.frame(x = rnorm(size), y = rnorm(size), z = rnorm(size)))
microbenchmark::microbenchmark(invisible(clustertree:::dt_knn(xyz, k = 30L, bkt_size = 30L, prune = FALSE)), times = 15L)
microbenchmark::microbenchmark(invisible(clustertree:::dt_knn(xyz, k = 30L, bkt_size = 30L, prune = TRUE)), times = 15L)
microbenchmark::microbenchmark(invisible(dbscan::kNN(xyz, k = 30L, sort = F, bucketSize = 30L)), times = 15L)

clustertree:::dt_knn(xyz, k = 30L, bkt_size = 30L, prune = TRUE)$dist[, -1] == dbscan::kNNdist(xyz, k = 30L, bucketSize = 30L)
invisible(clustertree:::dt_knn(xyz, k = 30L, bkt_size = 30L, prune = FALSE))

*/
