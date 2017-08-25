#ifndef DUAL_TREE_H
#define DUAL_TREE_H

#include <Rcpp.h>
using namespace Rcpp;

// ANN library + dual tree extensions
#include <ANN/ANN.h>
#include <ANN/kd_tree/kd_split.h> // ANNsplitRule
#include <DT/ANNdt/ANNkd_tree_dt.h>

// Utilities
#include "net_sort.h" // 2, 3, and 4 size sorting networks
#include <utilities.h> // R_INFO, profiling mode, etc.
#include <metrics.h> // different metrics

// STL or C lib includes
#include <algorithm> // std::for_each
#include <unordered_map> // unordered_map
#include <cassert> // assert
#include <queue> // queue

// Utility macros
// NOTE: Use caution when using, x is dereferenced in the processed, so expects a pointer
#include <typeinfo>
#define IS_SPLIT(x) (typeid(*x) == typeid(ANNkd_split))
#define AS_SPLIT(x) static_cast<ANNkd_split*>(x)
#define IS_LEAF(x) (typeid(*x) == typeid(ANNkd_leaf))
#define AS_LEAF(x) static_cast<ANNkd_leaf*>(x)

// Type definitions
// typedef std::pair<ANNkd_node*, ANNkd_node*> NODE_PAIR; // query node is always assumed as the first node

struct candidate_pair {
  bool beenChecked;
  ANNdist eps;
  candidate_pair() : beenChecked(false), eps(ANN_DIST_INF) {}
};

// To use as a default parameter
static NumericMatrix emptyMatrix = NumericMatrix(0, 0);
static List default_params = Rcpp::List::create(Rcpp::Named("bucketSize") = 10,
                                                Rcpp::Named("splitRule") = 5);

// ---- DualTree class definition ----
class DualTree {
public:
  Metric& m_dist;
  bool use_pruning; // whether to use a pruning-based approach
  int d; // dimension
  ANNkd_tree* qtree, *rtree; // query and reference tree pointers; could also be pointers to derived ANNkd_tree_dt types
  std::unordered_map<ANNkd_node*, const Bound& > bounds; // Various bounds per-node to fill in (node_ptr -> bounds)
  std::map< std::pair<int, int>, candidate_pair> BC_check; // Base Case point pair check: should default to false if no key found!

  #ifdef NDEBUG
    std::map<ANNkd_node*, char> node_labels; // in debug mode the nodes are labeled with characters
  #endif

  // Setup everything in the constructor!
  DualTree(const NumericMatrix& q_x, Metric& m, NumericMatrix& r_x = emptyMatrix, List& config = default_params);

  // Constructs one ANNkd_tree object. If pruning is set in the constructor, additional bounds are computed.
  virtual ANNkd_tree* ConstructTree(const NumericMatrix& x, const int bkt_sz = 15, ANNsplitRule = ANN_KD_SUGGEST);

  // Utility
  void PrintTree(ANNbool with_pts, bool);
  void printNode(ANNkd_split* N, int level);
  void printNode(ANNkd_leaf* N, int level);
  IntegerMatrix getBaseCases();

  // Distance calculation
  // By default, relies on a supplied metric and does incremental distance calculations. If any of the
  // distances exceed the given threshold parameters eps1 and eps2, the distance calculation is stopped,
  // and infinity is returned. Can be inherited and overridden by a derived class if more flexiblity is
  // needed.
  virtual ANNdist computeDistance(const int q_idx, const int r_idx,
                          ANNdist eps1 = ANN_DIST_INF,
                          ANNdist eps2 = ANN_DIST_INF);

  // Application-specific virtual methods: must be defined by the child!
  virtual void DFS(ANNkd_node* N_q, ANNkd_node* N_r) = 0; // regular full DFS
  virtual void pDFS(ANNkd_node* N_q, ANNkd_node* N_r) = 0; // pruning DFS variant
  virtual ANNdist Score(ANNkd_node* N_q, ANNkd_node* N_r) = 0;
  virtual ANNdist B(ANNkd_node* N_q) = 0;  // Bound function


  // Base case functions. Actual implementations to be derived, however
  // the simple base case can be used determine when to use a specific base
  // case based on if the trees are identical
  virtual ANNdist BaseCase(ANNkd_node* N_q, ANNkd_node* N_r);
  virtual ANNdist BaseCaseIdentity(ANNkd_node* N_q, ANNkd_node* N_r) = 0;
  virtual ANNdist BaseCaseNonIdentity(ANNkd_node* N_q, ANNkd_node* N_r) = 0;

  // With every dual tree traversal, pairs of points need to be scored. This method prevents equivalent
  // combinations from being scored multiple times.
  inline const bool hasBeenChecked(const ANNidx q_idx, const ANNidx r_idx){
    // inserts new entry (w/ value false) if the value didn't exist
    std::pair<int, int> key = std::minmax(q_idx, r_idx);
    const bool pair_visited = BC_check[key].beenChecked;
    if (!pair_visited) {
      R_INFO("== Calling base case for: q = " << q_idx << ", r = " << r_idx << "\n")
      BC_check[key].beenChecked = true; // Update knowledge known about the nodes
    }
    return pair_visited;
  }

  // Reset whether the base cases have been compared before
  void resetBaseCases(){
    for(std::map< std::pair<int, int>, candidate_pair>::iterator it = BC_check.begin(); it != BC_check.end(); ++it){
      BC_check[it->first].beenChecked = false;
    }
  }

  // Future work
  inline const void updateEps(const ANNidx q_idx, const ANNidx r_idx, ANNdist eps){
    BC_check[std::minmax(q_idx, r_idx)].eps = eps;
  }

  // The minimum distance between two nodes. If the convex subsets / rectangular subset they effectively
  // represent overlap, they'll have a negative distance, so return 0. Otherwise, the minimum distance
  // is defined as the distance between the centroids of each box, minus the maximum descendent distance
  // from centroid.
  inline ANNdist min_dist(ANNkd_node* N_q, ANNkd_node* N_r){// assume query and reference
    if (N_q == N_r) return 0; //equivalent nodes conatin the same points

    // Retrieve the bounds
    const Bound& nq_bound = bounds[N_q];
    const Bound& nr_bound = bounds[N_r];

    // Compute the minimum box distance
    // TODO: Move to other metrics, save result in cache
    ANNdist min_dist = m_dist(nq_bound.centroid, nr_bound.centroid) - nq_bound.lambda - nr_bound.lambda;

    // If it's negative, then the boxes overlap, and the minimum distance between any two points in either branch is 0.
    return(min_dist < 0 ? 0 : min_dist);
  }
};

#endif