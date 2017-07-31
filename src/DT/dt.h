#ifndef DUAL_TREE_H
#define DUAL_TREE_H

#include <Rcpp.h>
using namespace Rcpp;

// ANN library + dual tree extensions
#include <ANN/ANN.h>
#include "ANNdt/ANNkd_tree_dt.h"
#include "net_sort.h" // 2, 3, and 4 size sorting networks

// Utilities
#include <utilities.h> // R_INFO, profiling mode, etc.

// STL or C lib includes
#include <algorithm> // std::for_each
#include <unordered_map> // unordered_map
#include <cassert> // assert
#include <queue> // queue

// Utility macros
#define IS_SPLIT(x) (bool) (dynamic_cast<ANNkd_split*>(x) != NULL)
#define AS_SPLIT(x) dynamic_cast<ANNkd_split*>(x)
#define IS_LEAF(x) (bool) (dynamic_cast<ANNkd_leaf*>(x) != NULL)
#define AS_LEAF(x) dynamic_cast<ANNkd_leaf*>(x)

// Type definitions
// typedef std::pair<ANNkd_node*, ANNkd_node*> NODE_PAIR; // query node is always assumed as the first node

// ---- DualTree class definition ----
class DualTree {
  enum TREE_TYPE { QUERY_TREE, REF_TREE };
public:
  const bool use_pruning; // whether to use a pruning-based approach
  const int d; // dimension
  ANNkd_tree* qtree, *rtree; // query and reference tree pointers; could also be pointers to derived ANNkd_tree_dt types
  std::unordered_map<ANNkd_node*, const Bound& >* bounds; // Various bounds per-node to fill in (node_ptr -> bounds)
  std::map< std::pair<int, int>, bool>* BC_check; // Node pair base case check: should default to false if no key found!

  #ifdef NDEBUG
    std::unordered_map<ANNkd_node*, char> node_labels; // in debug mode the nodes are labeled with characters
  #endif

  DualTree(const bool prune, const int dim);
  virtual void setup(ANNkd_tree* kd_treeQ, ANNkd_tree* kd_treeR);
  virtual void setRefTree(ANNkd_tree* ref_tree) { rtree = ref_tree; };
  virtual void setQueryTree(ANNkd_tree* query_tree) { qtree = query_tree; };

  // Constructs one ANNkd_tree object. If pruning is set in the constructor, additional bounds are computed.
  virtual ANNkd_tree* ConstructTree(ANNpointArray x, const int nrow, const int ncol,
                                    const int bkt_sz = 15, ANNsplitRule = ANN_KD_SUGGEST);

  // Utility
  void PrintTree(ANNbool with_pts, bool);
  void printNode(ANNkd_split* N, int level);
  void printNode(ANNkd_leaf* N, int level);

  // Application-specific virtual methods: must be defined by the child!
  virtual void DFS(ANNkd_node* N_q, ANNkd_node* N_r) = 0; // regular full DFS
  virtual void pDFS(ANNkd_node* N_q, ANNkd_node* N_r) = 0; // pruning DFS variant
  virtual ANNdist BaseCase(ANNpoint p_q, ANNpoint p_r, const int q_idx, const int r_idx, ANNkd_node* N_q = NULL) = 0;
  virtual ANNdist Score(ANNkd_node* N_q, ANNkd_node* N_r) = 0;
  virtual ANNdist B(ANNkd_node* N_q) = 0;  // Bound function

  // Accessors to allow sub classes to access kd tree internals
  ANNkd_node* getRoot(TREE_TYPE type) { return type == QUERY_TREE ? qtree->root : rtree->root; };


  // With every dual tree traversal, pairs of points need to be scored. This method prevents equivalent
  // combinations from being scored multiple times.
  inline const bool hasBeenChecked(const ANNidx q_idx, const ANNidx r_idx){
    // inserts new entry (w/ value false) if the value didn't exist
    const bool pair_visited = (*BC_check)[std::minmax(q_idx, r_idx)];
    if (!pair_visited) {
      R_INFO("Calling base case for: q = " << q_idx << ", r = " << r_idx << ")\n")
      (*BC_check)[std::minmax(q_idx, r_idx)] = true; // Update knowledge known about the nodes
    }
    return pair_visited;
  }

  inline ANNdist min_dist(ANNkd_node* N_q, ANNkd_node* N_r){// assume query and reference
    if (N_q == N_r) return 0; //equivalent nodes conatin the same points

    // Assert they exist
    assert(bounds->find(N_q) != bounds->end() && bounds->find(N_q) != bounds->end());

    // Retrieve the bounds
    const Bound& nq_bound = (*bounds)[N_q];
    const Bound& nr_bound = (*bounds)[N_r];

    // Compute the minimum box distance
    ANNdist min_dist = annDist(d, nq_bound.centroid, nr_bound.centroid) - nq_bound.lambda - nr_bound.lambda;

    // If it's negative, then the boxes overlap, and the minimum distance between any two points in either branch is 0.
    return(min_dist < 0 ? 0 : min_dist);
  }
};

#endif