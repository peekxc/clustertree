#ifndef GREEDY_PRIMS_H
#define GREEDY_PRIMS_H

#include <Rcpp.h>
using namespace Rcpp;

// ANN includes
#include <ANN/ANNx.h>
#include <ANN/kd_tree/kd_tree.h>
#include <ANN/kd_tree/kd_util.h> // annBoxDist
#include <ANN/structures/pr_queue_k.h> // annBoxDist

// STL and other includes
#include <utilities.h>
#include <stack>
#include <unordered_map>

// Utility macros
// NOTE: Use caution when using, x is dereferenced in the processed, so expects a pointer
#include <typeinfo>
#define IS_SPLIT(x) (typeid(*x) == typeid(ANNkd_split))
#define AS_SPLIT(x) static_cast<ANNkd_split*>(x)
#define IS_LEAF(x) (typeid(*x) == typeid(ANNkd_leaf))
#define AS_LEAF(x) static_cast<ANNkd_leaf*>(x)

struct simple_edge {
  int to;
  ANNdist weight;
  simple_edge() : weight(ANN_DIST_INF), to(-1){};
  simple_edge(int id, ANNdist dist): weight(dist), to(id) { }
};

simple_edge findNearest(const int q_idx, ANNkd_tree* qtree, std::vector<int> ignore, double eps = 0);

void doSomething(ANNkd_leaf* node);
void doSomething(ANNkd_split* node);

// Stack-based variant of KNN
void knn_stack(const int q_idx, // the query point (index)
               const int k, 		// number of near neighbors to return
               ANNidxArray			nn_idx,			// nearest neighbor indices (returned)
               ANNdistArray		dd,				// the approximate nearest neighbor
               ANNkd_tree* qtree, // The kd tree used to make the points
               double eps = 0);

#endif