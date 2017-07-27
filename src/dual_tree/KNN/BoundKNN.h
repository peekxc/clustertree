#ifndef BOUND_KNN_H
#define BOUND_KNN_H

#include "../../ANN/ANNx.h"

// Pruning DT-DFS only: using regular Bound object for tree bounds, but
// provide additional bounds specific to KNN
struct BoundKNN {
  ANNdist B;                    // Actual computed bound
  unsigned int knn_known;       // the number of descendent (or held) points that have at least 1 non-inf KNN distance
  ANNdist min_knn;              // minimum distance to the kth nearest neighbor of the points within the node
  ANNdist max_real_knn;         // maximum non-inf knn distance
  BoundKNN() : B(ANN_DIST_INF), knn_known(0), min_knn(ANN_DIST_INF), max_real_knn(0) {};
  inline ANNdist maxKNN(const int npts){ return knn_known < npts ? ANN_DIST_INF : max_real_knn; }

  // Inserts point into knn priority queue.
  // inline void insert(ANNmin_k& qknn, ANNdist dist, ANNidx idx){
  //   qknn.insert(dist, idx);
  //   if (qknn.max_key() > max_real_knn){ max_real_knn = qknn.max_key(); } // Update max key if applicable
  // }
};

#endif