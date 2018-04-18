#ifndef DUALTREE_BORUVKA_H
#define DUALTREE_BORUVKA_H

#include "DualTreeSearch.h"
#include "simple_structs.h"
#include "union_find.h"
#include "profiler.h"

template <class METRIC_T>
class DualTreeBoruvka : public DualTreeSearch<METRIC_T> {
private:
  typedef DualTreeSearch<METRIC_T> dts; // dual tree search parent/base type definition
public:
  UnionFind CC; // Track connected components
  std::vector<double_edge> N; // map from component i -> j where j is the candidate nearest neighbor of i in a disjoint component
  std::vector<int> ALL_CC_SAME; // Are all descendent points of a given node in the same component?

  // Profiling
  int n_accesses;
  // std::unordered_map< pair<int,int>, ANNdist > dist_map;

  // Constructor
  DualTreeBoruvka(ANNkd_tree* qtree, ANNkd_tree* rtree, const NODE_INFO& qinfo, const NODE_INFO& rinfo, METRIC_T& m);

  // MST Base Case: Compute all the pairwise distances for each pair of points between the leaf the nodes
  void BaseCase(ANNkd_leaf* qn, ANNkd_leaf* rn);
  void BaseCaseCached(ANNkd_leaf* qn, ANNkd_leaf* rn);

  // MST Score Function
  double Score(ANNkd_node* qn, ANNkd_node* rn);

  // MST Update function for derived nodes
  void UpdateBounds(ANNkd_split* node, NODE_INFO& ninfo) override;
  void UpdateBounds(ANNkd_leaf* node, NODE_INFO& ninfo) override;

  NumericMatrix MST();

  // Convert the results in the priority queue suitable for returning over to R (1-based)
  NumericMatrix getResults();

  // Simple check to see if the components are fully connected, indicating a spanning tree is complete
  inline bool fully_connected(){
    bool is_fully_connected = true;
    for (int i = 0; i < CC.size - 1 && is_fully_connected; ++i){
      is_fully_connected = CC.Find(i) == CC.Find(i+1);
    }
    return is_fully_connected;
  }

  // Enable sorting of edges
  // static inline bool edge_comparator(const edge& e1, const edge& e2) { return e1.weight < e2.weight; }
};

#endif
