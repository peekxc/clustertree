#include <Rcpp.h>
using namespace Rcpp;

#include "dt_abstract.h"
#include "kd_tree.h"


DT_Abstract::DT_Abstract(ANNkd_tree* ref_tree, ANNkd_tree* query_tree) {
  rtree = ref_tree;
}

IntegerVector DT_Abstract::getIDXArray(){
  IntegerVector ids = Rcpp::no_init(rtree->n_pts);
  for (int i = 0; i < rtree->n_pts; ++i){
    ids[i] = rtree->pidx[i];
  }
  return ids;
}

IntegerVector DT_Abstract::child_nodes(){
  std::vector<int>* child_ids = new std::vector<int>();
  rtree->root->child_ids(*child_ids);

  // ANNkd_split* ref_root = (ANNkd_split*) rtree->root;
  //  ref_root_node = ref_root;
  // ref_root->
  return(wrap(*child_ids));
}

// std::vector<int>* DT_Abstract::child_nodes(ANNkd_node* node){
//   std::vector<int>* node_children = new std::vector<int>();
//   node->child_ids(*node_children); // adds by reference
//   return node_children;
// }

double DT_Abstract::min_dist(ANNkd_node* N_i, ANNkd_node* N_j){
  std::vector<int>* i_desc_ids = new std::vector<int>();
  std::vector<int>* j_desc_ids = new std::vector<int>();
  N_i->child_ids(*i_desc_ids);
  N_j->child_ids(*j_desc_ids);

  double min_dist = std::numeric_limits<ANNdist>::infinity(), dist;
  for (std::vector<int>::iterator i = i_desc_ids->begin(); i != i_desc_ids->end(); ++i){
    for (std::vector<int>::iterator j = j_desc_ids->begin(); j != j_desc_ids->end(); ++j){
      ANNpoint p_i = this->rtree->pts[*i];
      ANNpoint p_j = this->qtree->pts[*j];
      if ((dist = annDist(rtree->theDim(), p_i, p_j)) < min_dist){
        min_dist = dist;
      }
    }
  }
  return min_dist;
}

double DT_Abstract::max_dist(ANNkd_node* N_i, ANNkd_node* N_j){
  std::vector<int>* i_desc_ids = new std::vector<int>();
  std::vector<int>* j_desc_ids = new std::vector<int>();
  N_i->child_ids(*i_desc_ids);
  N_j->child_ids(*j_desc_ids);

  double max_dist = 0, dist;
  for (std::vector<int>::iterator i = i_desc_ids->begin(); i != i_desc_ids->end(); ++i){
    for (std::vector<int>::iterator j = j_desc_ids->begin(); j != j_desc_ids->end(); ++j){
      ANNpoint p_i = this->rtree->pts[*i];
      ANNpoint p_j = this->qtree->pts[*j];
      if ((dist = annDist(rtree->theDim(), p_i, p_j)) > max_dist){
        max_dist = dist;
      }
    }
  }
  return max_dist;
}

#include "kd_util.h"
NumericVector DT_Abstract::convex_subset(){
  NumericVector res = NumericVector(rtree->n_pts);
  for (int i = 0; i < rtree->n_pts; ++i){
    ANNdist box_dist = annBoxDistance(rtree->pts[i], rtree->bnd_box_lo,  rtree->bnd_box_lo, rtree->dim);
    res[i] = box_dist;
  }
  return res;
}

// double DT_Abstract::max_child_dist(ANNkd_node* N_i){
//
// }

