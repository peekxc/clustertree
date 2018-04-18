#include "RcppHeader.h"
#include "utilities.h"
#include "pr_queue_k.h"

#define min(x, y) (x<y?x:y)
#define max(x, y) (x<y?y:x)

// Given a 'dist' object and a integer k, return the k-th nearest neighbor for each point
// [[Rcpp::export]]
NumericVector knn_dist(const NumericVector& dist_x, const int k) {
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


#include <queue>

// [[Rcpp::export]]
NumericVector knn_dist2(const NumericVector& dist_x, const int k) {
  const int n = as<int>(dist_x.attr("Size"));

  // Create a set of priority queues
  std::vector< ANNmin_k* > knn_pq = std::vector<ANNmin_k*>(n);
  for (int i = 0; i < n; ++i){ knn_pq[i] = new ANNmin_k(k); }

  // Iterate through all of the distances
  int c = 0;
  for (NumericVector::const_iterator dist_ij = dist_x.begin(); dist_ij != dist_x.end(); ++dist_ij, ++c){
    const int i = INDEX_TO(c, n), j = INDEX_FROM(c, n, i);
    if ((*dist_ij) < knn_pq[i]->max_key()){ knn_pq[i]->insert(*dist_ij, j); }
    if ((*dist_ij) < knn_pq[j]->max_key()){ knn_pq[j]->insert(*dist_ij, i); }
  }

  // Get the max keys
  NumericVector knn_dist = Rcpp::no_init(n);
  for (int i = 0; i < n; ++i){ knn_dist[i] = knn_pq[i]->max_key(); }
  std::for_each(knn_pq.begin(), knn_pq.end(), [=](ANNmin_k* pq){ delete pq; }); // no longer need the priority queues

  // Return result
  return(knn_dist);
}


/*** R

x <- cbind(rnorm(10000), rnorm(10000))
iris_dist <- dist(x)
test <- clustertree:::knn_dist(iris_dist, 15L)
test2 <- clustertree:::knn_dist2(iris_dist, 15L)
truth <- dbscan::kNNdist(x, k = 15)[, 15]
truth == test == test2

microbenchmark::microbenchmark({ test <- clustertree:::knn_dist(iris_dist, 15L) }, times = 5L)
microbenchmark::microbenchmark({ test <- clustertree:::knn_dist2(iris_dist, 15L) }, times = 5L)
microbenchmark::microbenchmark({ dbscan::kNNdist(x, k = 15)[, 15] }, times = 5L)
*/
