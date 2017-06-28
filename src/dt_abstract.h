#ifndef DT_ABSTRACT_H
#define DT_ABSTRACT_H

#include <Rcpp.h>
using namespace Rcpp;

#include "ANN/ANN.h"

class DT_Abstract {
  ANNkd_tree* rtree, *qtree;
public:
  DT_Abstract(ANNkd_tree* ref_tree, ANNkd_tree* query_tree = NULL);
  IntegerVector getIDXArray();
  IntegerVector child_nodes();
  std::vector<int>* child_nodes(ANNkd_node* node);
  double min_dist(ANNkd_node* N_i, ANNkd_node* N_j);
  double max_dist(ANNkd_node* N_i, ANNkd_node* N_j);
  double max_child_dist(ANNkd_node* N_i);
  NumericVector convex_subset();
};

#endif