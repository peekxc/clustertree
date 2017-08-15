#ifndef DT_DTB_H
#define DT_DTB_H

#include <Rcpp.h>
using namespace Rcpp;

#include <ANN/ANN.h> // ANN defs
#include <metrics.h> // Distance metrics
#include <utilities.h> // R_info, profiling macros, etc.
#include <DT/dt.h> // Dual tree class definitions
#include <DT/structures/union_find.h> // Disjoint set data structure
#include <R_kdtree.h> // KD tree R-accessible functions

// Extend the DualTreeKNN class
#include <DT/KNN/dt_knn.h>

class DualTreeBoruvka : public DualTreeKNN {
protected:
  UnionFind CC;
  double_edge* EL; // edge list mapping component indices (starting with 1-n) to their closest
  std::unordered_map<ANNkd_node*, bool> ALL_CC_SAME; // are all descendent points of a given node in the same component?
  //edge* knn;  // Map between point index and kNN priority queue (k = 1 in this case)
  // std::unordered_map<unsigned int, double_edge> N; // Map from component index to shortest edge between components
  // std::unordered_map<unsigned int, ANNdist> D; // Map from component index to shortest edge distance
public:
  // Main constructors
  DualTreeBoruvka(const bool prune, const int dim, const int n, Metric* m = NULL); // default constructor
  void setup(ANNkd_tree* kd_treeQ, ANNkd_tree* kd_treeR);

  // New methods for the derived class
  NumericMatrix DTB(NumericMatrix& x);

  // Overridden Base case and score functions
  ANNdist BaseCaseIdentity(ANNkd_node* N_q, ANNkd_node* N_r) override;
  inline ANNdist Score(ANNkd_node* N_q, ANNkd_node* N_r) override;

  // Base-class functions
  //ANNdist B(ANNkd_node* N_q);
  //void pDFS(ANNkd_node* N_q, ANNkd_node* N_r);
  //void DFS(ANNkd_node* N_q, ANNkd_node* N_r){};
};

#endif