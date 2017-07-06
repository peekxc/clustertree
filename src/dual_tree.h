#ifndef DUAL_TREE_H
#define DUAL_TREE_H

#include <Rcpp.h>
using namespace Rcpp;

#include "ANN/ANN.h"
#include "pr_queue_k.h"
#include "kd_tree.h"

// C++11 includes
#include <unordered_map>

// Every node has a set number of bounds related to it that allow pruning of branch queries,
// each of which should only need to be computed once. The Bound struct stores these bounds,
// allowing the recursion to be memoized.
struct Bound {
  ANNdist B, rho, lambda, min_knn;
  ANNpoint centroid;
  ANNorthRect* bnd_box; // TODO: Have this point to the box computed in the tree construction
  Bound() : B(-1.0), rho(-1.0), lambda(-1.0), min_knn(ANN_DIST_INF), centroid(NULL) , bnd_box(NULL) { }
};

class DualTree {
  bool knn_identity;
public:
  const int d;
  ANNkd_tree* qtree, *rtree;
  ANNkd_node* N_q_par, *N_r_par;

  // Various bounds per-node to fill in (node_ptr -> bounds)
  std::unordered_map<ANNkd_node*, Bound>* bounds;

  // Map between point index and kNN priority queue
  std::unordered_map<ANNidx, ANNmin_k*>* knn;

  // Unordered map tracking whether base case has been called for a given pair of nodes
  std::map< std::pair<int, int>, bool>* BC_check; // should default to false if no key found!

  // Score and BaseCase functions
  double Score(ANNkd_node* N_q, ANNkd_node* N_r);
  void BaseCase(ANNpoint p_q, ANNpoint p_r, const int q_idx, const int r_idx,
                ANNkd_node* N_q, ANNkd_node* N_r);

  DualTree(ANNkd_tree* ref_tree, ANNkd_tree* query_tree);
  void knn_initialize(const int k);
  void test_cases(List&);
  void test_cases(List&, ANNkd_node*, int, bool);
  double basicDFS(ANNkd_node* node);
  void KNN(int k, NumericMatrix& dists, IntegerMatrix& ids);
  ANNdist B(ANNkd_node* N_q);
  void DFS(ANNkd_node* N_q, ANNkd_node* N_r);
  IntegerVector getIDXArray();
  IntegerVector child_ids(bool ref_tree = true);
  ANNdist min_dist(ANNkd_node* N_i, ANNkd_node* N_j);
  ANNdist max_dist(ANNkd_node* N_i, ANNkd_node* N_j);
  ANNdist max_child_dist(ANNkd_node* N_i, Bound& ni_bnd, bool ref_tree = true);
  ANNdist max_desc_dist(ANNkd_node* N_i, Bound& ni_bnd, bool ref_tree = true);
  ANNpoint centroid(ANNkd_node* N_i, Bound& ni_bnd, bool ref_tree = true);
  ANNorthRect& convex_subset(ANNkd_node* N_i, Bound& ni_bnd, bool ref_tree = true); // Keep return by reference
};


// Abstract base class
// struct ScoreFunction {
//   const DT_Abstract& env;
//   ScoreFunction(DT_Abstract dt_env) : env(dt_env) {};
//   virtual double operator() (ANNkd_node* N_q, ANNkd_node* N_r) = 0;
// };

#endif