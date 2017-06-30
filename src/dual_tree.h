#ifndef DUAL_TREE_H
#define DUAL_TREE_H

#include <Rcpp.h>
using namespace Rcpp;

#include "ANN/ANN.h"
#include "pr_queue_k.h"

// Forward Declaration
//struct ScoreFunction;

struct Bounds {
  std::map<ANNkd_node*, ANNdist > rho, lambda, d_min, d_max;
  // Bounds(){};
  // ANNdist operator()(ANNkd_node* node){
  //   if ( != mymap.end()))
  //   if (rho.at(node) > ){
  //
  //   }
  // }
};

struct Bound { ANNdist rho, lambda, d_min, d_max; };

class DualTree {
  ANNkd_tree* qtree, *rtree;
  ANNkd_node* N_q_par, *N_r_par;

  // Various bounds per-node to fill in (node_ptr -> bounds)
  std::map<ANNkd_node*, Bound>* bounds;

  // Map between point index and kNN priority queue
  std::map<ANNidx, ANNmin_k*>* knn;

  // Score and BaseCase functions
  double Score(ANNkd_node* N_q, ANNkd_node* N_r);
  double BaseCase(ANNpoint p_q, ANNpoint p_r, const int r_idx);
public:
  DualTree(ANNkd_tree* ref_tree, ANNkd_tree* query_tree, int k);
  void KNN(int k);
  double B1(ANNkd_node* N_q);
  double B2(ANNkd_node* N_q);
  void DFS(ANNkd_node* N_q, ANNkd_node* N_r);
  IntegerVector getIDXArray();
  IntegerVector child_ids(bool ref_tree = true);
  double min_dist(ANNkd_node* N_i, ANNkd_node* N_j);
  double max_dist(ANNkd_node* N_i, ANNkd_node* N_j);
  double max_child_dist(ANNkd_node* N_i);
  double max_desc_dist(ANNkd_node* N_i, bool ref_tree = true);
  NumericVector convex_subset(ANNkd_node* N_i, bool ref_tree = true);
};



// Abstract base class
// struct ScoreFunction {
//   const DT_Abstract& env;
//   ScoreFunction(DT_Abstract dt_env) : env(dt_env) {};
//   virtual double operator() (ANNkd_node* N_q, ANNkd_node* N_r) = 0;
// };
#endif