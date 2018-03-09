#include <Rcpp.h>
using namespace Rcpp;
#include "utilities.h"
#include "pr_queue_k.h"

#define min(x, y) (x<y?x:y)
#define max(x, y) (x<y?y:x)

// Given a 'dist' object and a integer k, return the k-th nearest neighbor for each point
// [[Rcpp::export]]
NumericVector knn_dist(NumericVector dist_x, const int k) {
  const int n = as<int>(dist_x.attr("Size"));
  NumericVector knn_dist = NumericVector(n);
  ANNmin_k* pr_queue;
  double cdist = 0.0;
  for (int i = 0; i < n; ++i){
    pr_queue = new ANNmin_k(k);
    for (int j = 0; j < n; ++j){
      if (i != j){
        cdist = dist_x[INDEX_TF(n, min(i, j), max(i, j))];
        pr_queue->insert(cdist, j);
      }
    }
    knn_dist[i] = pr_queue->max_key();
    delete pr_queue;
  }
  return(knn_dist);
}


/*** R

x <- cbind(rnorm(1000), rnorm(1000))
iris_dist <- dist(x)
test <- knn_dist(iris_dist, 15L)

truth <- apply(dbscan::kNNdist(x, k = 15), 1, max)
truth == test

*/
