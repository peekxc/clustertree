#include <Rcpp.h>
using namespace Rcpp;

#include "dual_tree.h"
#include "kd_tree.h"

// // static variables to update in the recursion
// #include "pr_queue_k.h"
// ANNmin_k* k_min_pts; // set of k closest points
//


// DualTree enables the dual-traversal of two kd-tree's at the same time, and
// for computing the associated bounds that come with them
DualTree::DualTree(ANNkd_tree* ref_tree, ANNkd_tree* query_tree, int k) {
  if (ref_tree->theDim() != query_tree->theDim()){
    stop("Dimensionality of the query set does not match the reference set.");
  }
  rtree = ref_tree, qtree = query_tree;
  N_q_par = N_r_par = NULL;
}

void DualTree::KNN(int k){

  // Initialize map of k-nearest neighbor queues
  knn = new std::map<ANNidx, ANNmin_k*>();
  for (int i = 0; i < qtree->n_pts; ++i){
   knn->at(qtree->pidx[i]) = new ANNmin_k(k);		// create set for closest k points
  }

  // Create map of bounds to save recursively
  bounds = new std::map<ANNkd_node*, Bound>();

  // start searching
  DFS(rtree->root, qtree->root);

  // deallocate closest point set
  // delete k_min_pts;

  return;
}

double DualTree::B1(ANNkd_node* N_q){
  std::vector<int> child_idx = std::vector<int>();
  N_q->child_ids(child_idx);
  double bound = 0, max_k;
  for (std::vector<int>::iterator idx = child_idx.begin(); idx != child_idx.end(); ++idx){
    max_k = knn->at(*idx)->max_key();
    bound = max_k > bound ? max_k : bound;
  }
  // Use infinity instead of max double
  bound = bound == ANN_DIST_INF ? std::numeric_limits<ANNdist>::infinity() : bound;
  return(bound);
}

// KNN score function
// struct ScoreKNN : ScoreFunction {
//     ScoreKNN(){};
//     double operator()(ANNkd_node* N_q, ANNkd_node* N_r){
//       N_q->max_child_dist()
//     }
// };


// KNN Score function
ANNdist DualTree::Score(ANNkd_node* N_q, ANNkd_node* N_r){
  // N_q->child_nodes()
  // ANNdist bound1 = std::max((ANNdist) k_min_pts->max_key(), );
  // if (min_dist(N_q, N_r) < bound){
  //   return (min_dist(N_q, N_r));
  // }
  return std::numeric_limits<ANNdist>::infinity();
}

// key field is distance
// info field is integer (idx)
ANNdist DualTree::BaseCase(ANNpoint p_q, ANNpoint p_r, const int r_idx){
  ANNdist d = annDist(this->rtree->theDim(), p_q, p_r);
  // if (d < k_min_pts->max_key())// TODO: make dynamic programming table for BASECASE pairs
  // { k_min_pts->insert(d, r_idx); }
  return d;
}

void DualTree::DFS(ANNkd_node* N_q, ANNkd_node* N_r){
  // ANNkd_split* split_query = dynamic_cast<ANNkd_split*>(N_q);
  // ANNkd_split* split_ref = dynamic_cast<ANNkd_split*>(N_r);

  if (Score(N_q, N_r) == std::numeric_limits<ANNdist>::infinity()){
    return;
  }

  // Get the points held in both the query and reference nodes
  std::vector<int>* query_pts = new std::vector<int>();
  std::vector<int>* ref_pts = new std::vector<int>();
  N_q->held_in_node(*query_pts), N_r->held_in_node(*ref_pts);
  for (std::vector<int>::iterator q_idx = query_pts->begin(); q_idx != query_pts->end(); ++q_idx){
    for (std::vector<int>::iterator r_idx = ref_pts->begin(); r_idx != ref_pts->end(); ++r_idx){
      BaseCase(qtree->pts[*q_idx], rtree->pts[*r_idx], *r_idx);
    }
  }

  // Get immediate children kd_nodes of the current query and references nodes
  std::vector<ANNkd_node*>* query_nodes = new std::vector<ANNkd_node*>();
  std::vector<ANNkd_node*>* ref_nodes = new std::vector<ANNkd_node*>();
  N_q->child_nodes(*query_nodes), N_r->child_nodes(*ref_nodes);

  // There are at minimum 0, and at maximum 2 nodes per node vector
  for (std::vector<ANNkd_node*>::iterator n_qc = query_nodes->begin(); n_qc != query_nodes->end(); ++n_qc){
    for (std::vector<ANNkd_node*>::iterator n_rc = ref_nodes->begin(); n_rc != ref_nodes->end(); ++n_rc){
      DFS((ANNkd_node*) *n_qc, (ANNkd_node*) *n_rc);
    }
  }
}

IntegerVector DualTree::getIDXArray(){
  IntegerVector ids = Rcpp::no_init(rtree->n_pts);
  for (int i = 0; i < rtree->n_pts; ++i){
    ids[i] = rtree->pidx[i];
  }
  return ids;
}

double DualTree::max_desc_dist(ANNkd_node* N_i, bool ref_tree){
  return 0;
}

// std::vector<int>* DualTree::child_nodes(ANNkd_node* node){
//   std::vector<int>* node_children = new std::vector<int>();
//   node->child_ids(*node_children); // adds by reference
//   return node_children;
// }

double DualTree::min_dist(ANNkd_node* N_i, ANNkd_node* N_j){

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

double DualTree::max_dist(ANNkd_node* N_i, ANNkd_node* N_j){
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

NumericVector DualTree::convex_subset(ANNkd_node* N_i, bool ref_tree){

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
  NumericVector res = NumericVector::create(*bnd_box.lo, *bnd_box.hi);
  return res;
}