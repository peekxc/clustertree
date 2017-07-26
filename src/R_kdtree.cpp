#include <Rcpp.h>
using namespace Rcpp;

// ANN library
#include "ANN/ANN.h"
#include "kd_pr_search.h"

// [[Rcpp::export]]
List kd_knn(NumericMatrix query_x, SEXP tree_ptr, int k, bool priority){

  // Extract the kdtree back
  Rcpp::XPtr<ANNkd_tree> kd_tree_ptr(tree_ptr);
  ANNkd_tree& kd_tree = *kd_tree_ptr;

  // Note: the search also returns the point itself (as the first hit)!
  // So we have to look for k+1 points.
  ANNdistArray dists = new ANNdist[k+1];
  ANNidxArray nnIdx = new ANNidx[k+1];

  // Distance matrix of kNN distances
  NumericMatrix d(query_x.nrow(), k);

  // Id matrix of knn indices
  IntegerMatrix id(query_x.nrow(), k);

  // Loop through the query points and perform the search
  for (int i=0; i < query_x.nrow(); i++) {
    if (!(i % 100)) Rcpp::checkUserInterrupt();

    // Retrieve the query point
    NumericVector query_pt = query_x.row(i);
    ANNpoint queryPt = &query_pt.at(0);


    #ifdef ANN_PERF
        annResetCounts();			// reset counters
    #endif

    // Search the kd tree
    kd_tree.annkSearch(queryPt, k+1, nnIdx, dists, 0);

    #ifdef ANN_PERF
      annUpdateStats();
    #endif
    // kd_tree.annkPriSearch();

    // Remove self matches
    IntegerVector ids = IntegerVector(nnIdx, nnIdx+k+1);
    LogicalVector take = ids != i;
    ids = ids[take];
    id(i, _) = ids + 1;

    // Convert to regular, non-squared euclidean distances
    NumericVector ndists = NumericVector(dists, dists+k+1)[take];
    d(i, _) = sqrt(ndists);
  }
  List res = List::create(_["dist"] = d, _["id"] = id);
  return(res);
}

// [[Rcpp::export]]
List kdtree(NumericMatrix x, const int bkt_size) {

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

  // Cleanup
  // delete kdTree;
  // delete [] dists;
  // delete [] nnIdx;
  // annDeallocPts(dataPts);
  annClose();

  // prepare results
  Rcpp::XPtr<ANNkd_tree> p(kdTree, true);
  List ret = List::create(_["kdtree_ptr"] = p);
  return ret;
}