#ifndef DUAL_TREE_H
#define DUAL_TREE_H

#include <Rcpp.h>
using namespace Rcpp;

#include "ANN/ANN.h"
#include "pr_queue_k.h"
#include "kd_tree.h"
#include "kd_util.h"

#include <algorithm> // std::for_each
#include <unordered_map> // unordered_map
#include <cassert> // assert
#include <queue> // queue

//#define NDEBUG 1 // <-- for 'debug' mode, will print out info to R session
#undef NDEBUG // <-- for 'production' mode, will remove IO

#ifdef NDEBUG
  #define ANN_PERF // Compile ANN library for performance testing
  #define R_INFO(x) Rcpp::Rcout << x;
  #define R_PRINTF(x, ...) Rprintf(x, __VA_ARGS__);
#else
  #define R_INFO(x)
  #define R_PRINTF(x, ...)
#endif

// Utility macros
#define IS_SPLIT(x) (bool) (dynamic_cast<ANNkd_split*>(x) != NULL)
#define AS_SPLIT(x) dynamic_cast<ANNkd_split*>(x)
#define IS_LEAF(x) (bool) (dynamic_cast<ANNkd_leaf*>(x) != NULL)
#define AS_LEAF(x) dynamic_cast<ANNkd_leaf*>(x)

#define VI(x) std::vector<x>::iterator

// Type definitions
typedef std::pair<ANNkd_node*, ANNkd_node*> NODE_PAIR; // query node is always assumed as the first node

// Every node has a set number of bounds related to it that allow pruning of branch queries,
// each of which should only need to be computed once. The Bound struct stores these bounds,
// allowing the recursion to be memoized.
struct Bound {
  ANNdist B, rho, lambda;
  ANNpoint centroid;
  ANNorthRect* bnd_box; // TODO: Have this point to the box computed in the tree construction
  Bound() : B(-1.0), rho(-1.0), lambda(-1.0), centroid(NULL) , bnd_box(NULL) { }
};

// ---- DualTree class definition ----
class DualTree {
protected:
  const int d; // dimension
  ANNkd_tree* qtree, *rtree; // query and reference tree pinters
  std::map< std::pair<int, int>, bool>* BC_check; // Node pair base case check: should default to false if no key found!
public:
  // Constructors and destructors: just need two trees
  DualTree(ANNkd_tree* ref_tree, ANNkd_tree* query_tree);
  // ~DualTree(); // virtual destructor

  // Utility
  void PrintTree(bool);
  // void PrintStatistics(bool);
  IntegerVector getIDXArray();
  IntegerVector child_ids(bool ref_tree = true);

  // Application-specific virtual methods: must be defined by the child!
  virtual void DFS(ANNkd_node* N_q, ANNkd_node* N_r) = 0; // regular full DFS
  virtual void pDFS(ANNkd_node* N_q, ANNkd_node* N_r) = 0; // pruning DFS variant
  virtual void BaseCase(ANNpoint p_q, ANNpoint p_r, const int q_idx, const int r_idx, ANNkd_node* N_q = NULL) = 0;
  virtual ANNdist Score(ANNkd_node* N_q, ANNkd_node* N_r) = 0;
  virtual ANNdist B(ANNkd_node* N_q) = 0;  // Bound function

  // Tree-specific virtual functions: while they can be re-derived, polymorphism is handled internally
  // through function overloading via the ANNkd_node class and it's derivatives, so in many cases the default
  // functionality provided by this class is suitable
  virtual ANNdist min_dist(ANNkd_node* N_i, ANNkd_node* N_j);
  virtual ANNdist max_dist(ANNkd_node* N_i, ANNkd_node* N_j);
  virtual ANNdist max_child_dist(ANNkd_node* N_i, Bound& ni_bnd, bool ref_tree = true);
  virtual ANNdist max_desc_dist(ANNkd_node* N_i, Bound& ni_bnd, bool ref_tree = true);
  virtual ANNpoint centroid(ANNkd_node* N_i, Bound& ni_bnd, bool ref_tree = true);
  virtual ANNorthRect& convex_subset(ANNkd_node* N_i, Bound& ni_bnd, bool ref_tree = true); // Keep ANNorthRect return by reference

  // void test_cases(List&);
  // void test_cases(List&, ANNkd_node*, int, bool);
};

#endif