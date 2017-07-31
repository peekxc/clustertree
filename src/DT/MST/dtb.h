#ifndef DT_DTB_H
#define DT_DTB_H

#include <Rcpp.h>
using namespace Rcpp;

#include <ANN/ANN.h> // ANN defs
#include <metrics.h> // Distance metrics
#include <utilities.h> // R_info, profiling macros, etc.
#include "../dt.h" // Dual tree class definitions
#include <DT/structures/union_find.h> // Disjoint set data structure
#include <R_kdtree.h> // KD tree R-accessible functions

// Borrow some tools from KNN
#include <DT/KNN/BoundKNN.h>

class DualTreeBoruvka : public DualTree {
protected:
  ANNkd_node* N_q_par, *N_r_par;
  std::unordered_map<unsigned int, double_edge> N; // Map from component index to shortest edge between components
  std::unordered_map<unsigned int, ANNdist> D; // Map from component index to shortest edge distance
  UnionFind CC;

  // Tools from DT KNN
  std::unordered_map<ANNidx, edge>* knn;  // Map between point index and kNN priority queue (k = 1 in this case)
  std::unordered_map<ANNkd_node*, BoundKNN& >* bnd_knn; // Useful node-wide KNN-specific bounds/properties
public:
  // Main constructors
  DualTreeBoruvka(const bool prune, const int dim, const int n); // default constructor
  virtual void setup(ANNkd_tree* kd_treeQ, ANNkd_tree* kd_treeR);

  // New methods for the derived class
  void DTB(NumericMatrix& x);

  // Base-class functions
  virtual ANNdist B(ANNkd_node* N_q);
  virtual void pDFS(ANNkd_node* N_q, ANNkd_node* N_r){};
  virtual void DFS(ANNkd_node* N_q, ANNkd_node* N_r){};

  // Score and BaseCase functions: try to inline if possible
  virtual inline ANNdist Score(ANNkd_node* N_q, ANNkd_node* N_r);
  virtual inline ANNdist BaseCase(ANNpoint p_q, ANNpoint p_r, const int q_idx, const int r_idx, ANNkd_node* N_q);
};

#endif