// R_kNN.cpp
// R-facing code to perform k-nearest neighbor.
// Single tree kNN function copied with permission from Michael Hahsler's code in the 'dbscan' package (GPLv3).
#include <Rcpp.h>
using namespace Rcpp;

#include <utilities.h> // R_INFO, profiling mode, etc. must include first
#include <DT/KNN/dt_knn.h> // Dual Tree KNN implementation
#include "R_kdtree.h" // R Interface to ANN KD trees
#include <ANN/ANN_util.h>
#include "R_kNN.h"

using namespace Rcpp;

// returns knn + dist
List kNN_int(NumericMatrix data, int k,
  int type, int bucketSize, int splitRule, double approx) {

  // copy data
  int nrow = data.nrow();
  int ncol = data.ncol();
  ANNpointArray dataPts = annAllocPts(nrow, ncol);
  for(int i = 0; i < nrow; i++){
    for(int j = 0; j < ncol; j++){
      (dataPts[i])[j] = data(i, j);
    }
  }
  //Rprintf("Points copied.\n");

  // create kd-tree (1) or linear search structure (2)
  ANNpointSet* kdTree = NULL;
  if (type==1){
    kdTree = new ANNkd_tree(dataPts, nrow, ncol, bucketSize,
      (ANNsplitRule)  splitRule);
  } else{
    kdTree = new ANNbruteForce(dataPts, nrow, ncol);
  }
  //Rprintf("kd-tree ready. starting DBSCAN.\n");

  NumericMatrix d(nrow, k);
  IntegerMatrix id(nrow, k);

  // Note: the search also returns the point itself (as the first hit)!
  // So we have to look for k+1 points.
  ANNdistArray dists = new ANNdist[k+1];
  ANNidxArray nnIdx = new ANNidx[k+1];

  for (int i=0; i<nrow; i++) {
    if (!(i % 100)) Rcpp::checkUserInterrupt();

    ANNpoint queryPt = dataPts[i];

    if(type==1) kdTree->annkSearch(queryPt, k+1, nnIdx, dists, approx);
    else kdTree->annkSearch(queryPt, k+1, nnIdx, dists);

    // remove self match
    IntegerVector ids = IntegerVector(nnIdx, nnIdx+k+1);
    LogicalVector take = ids != i;
    ids = ids[take];
    id(i, _) = ids + 1;

    NumericVector ndists = NumericVector(dists, dists+k+1)[take];
    d(i, _) = sqrt(ndists);
  }

  // cleanup
  delete kdTree;
  delete [] dists;
  delete [] nnIdx;
  annDeallocPts(dataPts);
  annClose();

  // prepare results
  List ret;
  ret["dist"] = d;
  ret["id"] = id;
  ret["k"] = k;
  return ret;
}

// R-facing API for the KNN-focused dual tree traversal
// TODO: approx unused right now
// [[Rcpp::export]]
List dt_kNN_int(const NumericMatrix& q_x, const int k, int bucketSize, int splitRule, SEXP metric_ptr, // The user-specified metric to use
                NumericMatrix r_x = NumericMatrix())
  {

  // Construct the dual tree KNN instance
  Metric& metric = getMetric(metric_ptr);
  List config = List::create(_["bucketSize"] = bucketSize, _["splitRule"] = splitRule);
  DualTreeKNN dt_knn = DualTreeKNN(q_x, metric, r_x, config);

  // Note: the search also returns the point itself (as the first hit)!
  // So we have to look for k+1 points.
  List res = dt_knn.KNN(k+1);

  // Return results
  return(res);
}

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



