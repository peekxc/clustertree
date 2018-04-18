// R_kNN.cpp
// R-facing code to perform k-nearest neighbor.
// Single tree kNN function copied and slightly editted with permission from Michael Hahsler's code in the 'dbscan' package (GPLv3).
#include "RcppHeader.h"
#include "utilities.h" // R_INFO, profiling mode, etc. must include first
#include "ANN_util.h" // matrixToANNpointArray
#include "R_kNN.h" // main header
#include "R_kdtree.h" // R Interface to ANN KD trees
// #include <DT/KNN/dt_knn.h> // Dual Tree KNN implementation

using namespace Rcpp;

// returns knn + dist
List kNN_int(const NumericMatrix& x, int k,
             int type, int bucketSize, int splitRule, double approx,
             NumericMatrix& r_x = emptyMatrix) {

  // Copy data over to ANN point array
  bool qr_same = true;
  ANNpointArray queryPts = matrixToANNpointArray(x);
  ANNpointArray refPts;
  if (r_x.size() <= 1){
    r_x = x;
    refPts = queryPts;
  } else {
    qr_same = false;
    refPts = matrixToANNpointArray(r_x);
  }

  // Query data dimensions
  int q_n = x.nrow(), q_d = x.ncol();
  int r_n = r_x.nrow(), r_d = r_x.ncol();
  if (q_d != r_d){ Rcpp::stop("Query and reference sets must have the same dimension."); }

  // Create the ANN kd tree
  ANNkd_tree* kdTree = new ANNkd_tree(refPts, r_n, r_d, bucketSize, (ANNsplitRule) splitRule);

  // Reserve the output
  NumericMatrix dist(q_n, k);
  IntegerMatrix id(q_n, k);

  // Note: the search also returns the point itself (as the first hit)!
  // So we have to look for k+1 points.
  ANNdistArray dists = new ANNdist[k+1];
  ANNidxArray nnIdx = new ANNidx[k+1];

  for (int i=0; i < q_n; i++) {
    if (!(i % 100)) Rcpp::checkUserInterrupt();

    ANNpoint queryPt = queryPts[i];
    kdTree->annkSearch(queryPt, k+1, nnIdx, dists, approx);

    // remove self match
    IntegerVector ids = IntegerVector(nnIdx, nnIdx+k+1);
    LogicalVector take = ids != i;
    ids = ids[take];
    id(i, _) = ids + 1;

    NumericVector ndists = NumericVector(dists, dists+k+1);
    dist(i, _) = sqrt(ndists[take]);
  }

  // cleanup
  delete kdTree;
  delete [] dists;
  delete [] nnIdx;
  annDeallocPts(queryPts);
  if (!qr_same) annDeallocPts(refPts);
  annClose();

  // prepare results
  List ret;
  ret["dist"] = dist;
  ret["id"] = id;
  ret["k"] = k;
  return ret;
}

// R-facing API for the KNN-focused dual tree traversal
// TODO: approx unused right now
// List dt_kNN_int(const NumericMatrix& q_x, const int k, int bucketSize, int splitRule, SEXP metric_ptr, // The user-specified metric to use
//                 NumericMatrix r_x = NumericMatrix())
//   {
//
//   // Construct the dual tree KNN instance
//   Metric& metric = getMetric(metric_ptr);
//   List config = List::create(_["bucketSize"] = bucketSize, _["splitRule"] = splitRule);
//   DualTreeKNN dt_knn = DualTreeKNN(q_x, metric, r_x, config);
//
//   // Note: the search also returns the point itself (as the first hit)!
//   // So we have to look for k+1 points.
//   List res = dt_knn.KNN(k+1);
//
//   // Return results
//   return(res);
// }

// Performance profiling
// void perf_test(){
//   List res;
// #ifdef PROFILING
//   annResetStats(r_x.nrow());
//
//   // Regular kdtree search
//   List ann_kdtree = kdtree(r_x, bucketSize); // create regular kdtree with the same bucket size
//   BEGIN_PROFILE()
//     List res1 = kd_knn(q_x, (SEXP) ann_kdtree["kdtree_ptr"], k, false);
//   REPORT_TIME("Regular KNN search")
//
//     // Print statistics
//     Rcout << "Regular KDtree search performance: " << std::endl;
//   annPrintStats((ANNbool) false);
//
//   // Dual tree equivalent
//   Rcout << "DualTree search performance: " << std::endl;
//   annResetStats(r_x.nrow());
//   annResetCounts();	// reset stats for a set of queries
//   BEGIN_PROFILE()
// }



