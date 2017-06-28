#ifndef DT_ABSTRACT_H
#define DT_ABSTRACT_H

#include <Rcpp.h>
using namespace Rcpp;

#include "ANN/ANN.h"

class DT_Abstract {
  ANNkd_tree* rtree, *qtree;
  double Score(ANNkd_node* N_q, ANNkd_node* N_r);
  double BaseCase(ANNpoint p_q, ANNpoint p_r);
public:
  DT_Abstract(ANNkd_tree* ref_tree, ANNkd_tree* query_tree = NULL);
  void DFS(ANNkd_node* N_q, ANNkd_node* N_r);
  IntegerVector getIDXArray();
  IntegerVector child_ids(bool ref_tree = true);
  std::vector<int>* child_nodes(ANNkd_node* node);
  double min_dist(ANNkd_node* N_i, ANNkd_node* N_j);
  double max_dist(ANNkd_node* N_i, ANNkd_node* N_j);
  double max_child_dist(ANNkd_node* N_i);
  double max_desc_dist(ANNkd_node* N_i, bool ref_tree = true);
  NumericVector convex_subset(ANNkd_node* N_i, bool ref_tree = true);
};

#endif