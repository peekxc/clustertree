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

struct edgeFT {
  unsigned int from, to;
  edgeFT() { from = -1; to = -1; } // default constructor
  edgeFT(int from_id, int to_id) : from(from_id), to(to_id) { }
};

class DualTreeBoruvka : public DualTreeKNN {
protected:
  UnionFind CC;
  std::vector<edgeFT> N; // component index map to edge of point in the component to the candidate nearest neighbor outside of the component
  std::vector<ANNdist> D; // Component index map to the distance between the two points above
  double_edge* EL; // edge list mapping component indices (starting with 1-n) to their closest

  // Are all descendent points of a given node in the same component?
  // -1 if no, o.w. stores the component index
  std::unordered_map<ANNkd_node*, int> ALL_CC_SAME;

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

  // Simple check to see if the components are fully connected, indicating a spanning tree is complete
  bool fully_connected(){
    bool is_fully_connected = true;
    for (int i = 0; i < CC.size - 1 && is_fully_connected; ++i){
      is_fully_connected = CC.Find(i) == CC.Find(i+1);
    }
    return is_fully_connected;
  }
};

#endif