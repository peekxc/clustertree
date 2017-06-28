#include <Rcpp.h>
using namespace Rcpp;

#include "dt_abstract.h"
#include "kd_tree.h"

// DT_Abstract does not act as an inheritance-based, object-oriented type class, but rather
// as an abstract class for enabling the dual-traversal of two kd-tree's at the same time, and
// for computing the associated bounds that come with them
DT_Abstract::DT_Abstract(ANNkd_tree* ref_tree, ANNkd_tree* query_tree) {
  rtree = ref_tree, qtree = query_tree;
}

double DT_Abstract::Score(ANNkd_node* N_q, ANNkd_node* N_r){
  return 0;
}
double DT_Abstract::BaseCase(ANNpoint p_q, ANNpoint p_r){
  return 0;
}

void DT_Abstract::DFS(ANNkd_node* N_q, ANNkd_node* N_r){
  ANNkd_split* split_query = dynamic_cast<ANNkd_split*>(N_q);
  ANNkd_split* split_ref = dynamic_cast<ANNkd_split*>(N_r);


  if (Score(N_q, N_r) == std::numeric_limits<ANNdist>::infinity()){
    return;
  }

  // Get the points held in both the query and reference nodes
  std::vector<int>* query_pts = new std::vector<int>();
  std::vector<int>* ref_pts = new std::vector<int>();
  N_q->held_in_node(*query_pts), N_r->held_in_node(*ref_pts);

  std::vector<int>* query_idx = new std::vector<int>();
  std::vector<int>* ref_idx = new std::vector<int>();
  N_q->child_ids(*query_pts), N_r->child_ids(*ref_pts);
  for (std::vector<int>::iterator n_qc = query_pts->begin(); n_qc != query_pts->end(); ++n_qc){
    for (std::vector<int>::iterator n_rc = ref_idx->begin(); n_rc != ref_idx->end(); ++n_rc){
      // DFS((ANNkd_node*) *n_qc, (ANNkd_node*) *n_rc);
    }
  }

}

IntegerVector DT_Abstract::getIDXArray(){
  IntegerVector ids = Rcpp::no_init(rtree->n_pts);
  for (int i = 0; i < rtree->n_pts; ++i){
    ids[i] = rtree->pidx[i];
  }
  return ids;
}

double DT_Abstract::max_desc_dist(ANNkd_node* N_i, bool ref_tree){
  return 0;
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
#include "ANN/ANNx.h"
#include "ANN/ANNperf.h"

NumericVector DT_Abstract::convex_subset(ANNkd_node* N_i, bool ref_tree){

  // Choose the tree
  ANNkd_tree* ctree = ref_tree ? rtree : qtree;

  // Get child ids
  std::vector<int>* n_child_ids = new std::vector<int>();
  N_i->child_ids(*n_child_ids);
  ANNidxArray pidx = (ANNidxArray) &n_child_ids[0];

  // Get the points associated with the child ids
  ANNpointArray pa[n_child_ids->size()];
  for (int i = 0; i < n_child_ids->size(); ++i){ pa[i] = (pa[pidx[(i)]]); }

  // Compute the bounding box
  ANNorthRect bnd_box(rtree->theDim());
  annEnclRect((ANNpointArray) pa,
              (ANNidxArray) pidx,
              (int) n_child_ids->size(),
              (int) rtree->theDim(),
              bnd_box); // construct bounding rectangle

  // Return results
  NumericVector res = NumericVector::create(bnd_box.lo, bnd_box.hi);
  return res;
}