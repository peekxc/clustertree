#include <Rcpp.h>
using namespace Rcpp;

#include "dual_tree.h"
#include "kd_tree.h"

// // static variables to update in the recursion
// #include "pr_queue_k.h"
// ANNmin_k* k_min_pts; // set of k closest points

// DualTree enables the dual-traversal of two kd-tree's at the same time, and
// for computing the associated bounds that come with them
DualTree::DualTree(ANNkd_tree* ref_tree, ANNkd_tree* query_tree, int k) : d(ref_tree->theDim()) {
  if (ref_tree->theDim() != query_tree->theDim()){
    stop("Dimensionality of the query set does not match the reference set.");
  }
  rtree = ref_tree, qtree = query_tree;
  N_q_par = N_r_par = NULL;
}

List DualTree::test_cases() {
  Rcout << "Beginning test cases" << std::endl;
  R_FlushConsole();

  List test_results = List();

  // Test getting the convex subset
  ANNorthRect r_bnd_box = this->convex_subset(rtree->root, true);
  ANNorthRect q_bnd_box = this->convex_subset(qtree->root, false);

  // Test retrieving the centroid
  NumericVector r_res = NumericVector(d), q_res = NumericVector(d);
  ANNpoint r_centroid = centroid(rtree->root, true), q_centroid = centroid(qtree->root, false);
  std::copy(r_centroid, r_centroid + d, r_res.begin());
  std::copy(q_centroid, q_centroid + d, q_res.begin());
  test_results["r_centroid"] = r_res;
  test_results["q_centroid"] = q_res;

  // Return results
  return test_results;
}

void DualTree::KNN(int k, NumericMatrix& dists, IntegerMatrix& ids){

  // Initialize map of k-nearest neighbor queues
  // TODO: Map back to integer vector
  knn = new std::map<ANNidx, ANNmin_k*>();
  for (int i = 0; i < qtree->n_pts; ++i) {
    // Create set for closest k points
    knn->insert(std::pair<ANNidx, ANNmin_k*>(qtree->pidx[i], new ANNmin_k(k)));
  }

  // Create map of bounds to memoize the results of the recursion
  bounds = new std::map<ANNkd_node*, Bound>();

  // Start the search!
  Rcout << "Starting DFS search!" << std::endl;
  DFS(rtree->root, qtree->root);

  // Copy over the distances and ids
  Rcout << "Copying distance and ids over to R memory" << std::endl;
  for (int i = 0; i < qtree->n_pts; ++i){
    for (int j = 0; j < k; j++) {		// extract the k-th closest points
      dists(i, j) = knn->at(qtree->pidx[i])->ith_smallest_key(i);
      ids(i, j) = knn->at(qtree->pidx[i])->ith_smallest_info(i);;
    }
  }

  // deallocate closest point set
  // delete k_min_pts;

  return;
}

double DualTree::B1(ANNkd_node* N_q){
  std::vector<int> child_idx = std::vector<int>();
  N_q->node_ids(child_idx);
  double bound = 0, max_k;
  for (std::vector<int>::iterator idx = child_idx.begin(); idx != child_idx.end(); ++idx){
    max_k = knn->at(*idx)->max_key();
    bound = max_k > bound ? max_k : bound;
  }

  // Get child nodes
  std::vector<ANNkd_node*> child_nodes = std::vector<ANNkd_node*>();
  N_q->child_nodes(child_nodes);
  for (std::vector<ANNkd_node*>::iterator node = child_nodes.begin(); node != child_nodes.end(); ++node){
    bound = std::max(bound, B1(*node));
  }

  // Use infinity instead of max double
  bound = bound == ANN_DIST_INF ? std::numeric_limits<ANNdist>::infinity() : bound;

  return(bound);
}

double DualTree::B2(ANNkd_node* N_q){
  std::vector<int> child_idx = std::vector<int>();
  N_q->node_ids(child_idx);
  double bound = std::numeric_limits<ANNdist>::infinity(), min_k;
  for (std::vector<int>::iterator idx = child_idx.begin(); idx != child_idx.end(); ++idx){
    min_k = knn->at(*idx)->ith_smallest_key(0);
    bound = min_k < bound ? min_k : bound;
  }
  bound += 2 * max_desc_dist(N_q, false);
  // Use infinity instead of max double
  bound = bound == ANN_DIST_INF ? std::numeric_limits<ANNdist>::infinity() : bound;

  return(bound);
}

// Recursive bound on a given query node
ANNdist DualTree::B(ANNkd_node* N_q){
  Rcout << "Computing the bound of a query node" << std::endl;

  ANNdist inf_const = std::numeric_limits<ANNdist>::infinity();
  ANNdist b1 = 0, b2 = inf_const, b3 = inf_const, b4 = inf_const;
  b1 = b2 = b3 = b4 = inf_const;

  // ----- First bound -----
  std::vector<int> query_ids = std::vector<int>();
  N_q->node_ids(query_ids);

  for (std::vector<int>::iterator idx = query_ids.begin(); idx != query_ids.end(); ++idx){
    ANNdist max_k = knn->at(*idx)->max_key();
    b1 = max_k > b1 ? max_k : b1;
  }

  std::vector<ANNkd_node*> child_nodes = std::vector<ANNkd_node*>();
  N_q->child_nodes(child_nodes);
  for (std::vector<ANNkd_node*>::iterator idx = child_nodes.begin(); idx != child_nodes.end(); ++idx){
    ANNdist child_bound = B(*idx);
    b1 = child_bound > b1 ? child_bound : b1;
  }

  // ----- Second bound -----
  for (std::vector<int>::iterator idx = query_ids.begin(); idx != query_ids.end(); ++idx){
    ANNdist child_bound = knn->at(*idx)->max_key() + max_child_dist(N_q, false) + max_desc_dist(N_q, false);
    b2 = child_bound < b2 ? child_bound : b2;
  }

  // ----- Third bound -----
  for (std::vector<ANNkd_node*>::iterator N_c = child_nodes.begin(); N_c != child_nodes.end(); ++N_c){
    ANNdist child_bound = B(*N_c) + 2 * (max_desc_dist(N_q) - max_desc_dist(*N_c));
    b3 = child_bound < b3 ? child_bound : b3;
  }

  // ----- Fourth bound -----
  b4 = N_q_par == NULL ? inf_const : B(N_q_par);

  // Final bound
  return std::min(std::min(b1, b2), std::min(b3, b4));
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
  if (min_dist(N_q, N_r) < B(N_q)){
    return (min_dist(N_q, N_r));
  }
  return std::numeric_limits<ANNdist>::infinity();
}

// key field is distance
// info field is integer (idx)
ANNdist DualTree::BaseCase(ANNpoint p_q, ANNpoint p_r, const int q_idx, const int r_idx){
  ANNdist dist = annDist(d, p_q, p_r);
  if (dist < knn->at(q_idx)->max_key())// TODO: make dynamic programming table for BASECASE pairs
  { knn->at(q_idx)->insert(d, r_idx); }
  return d;
}

void DualTree::DFS(ANNkd_node* N_q, ANNkd_node* N_r){
  // ANNkd_split* split_query = dynamic_cast<ANNkd_split*>(N_q);
  // ANNkd_split* split_ref = dynamic_cast<ANNkd_split*>(N_r);

  Rcout << "Computing Score" << std::endl;
  R_FlushConsole();
  if (Score(N_q, N_r) == std::numeric_limits<ANNdist>::infinity()){
    return;
  }

  // Get the points held in both the query and reference nodes
  Rcout << "Traversing Base cases" << std::endl;
  std::vector<int>* query_pts = new std::vector<int>();
  std::vector<int>* ref_pts = new std::vector<int>();
  N_q->node_ids(*query_pts), N_r->node_ids(*ref_pts);
  for (std::vector<int>::iterator q_idx = query_pts->begin(); q_idx != query_pts->end(); ++q_idx){
    for (std::vector<int>::iterator r_idx = ref_pts->begin(); r_idx != ref_pts->end(); ++r_idx){
      BaseCase(qtree->pts[*q_idx], rtree->pts[*r_idx], *q_idx, *r_idx);
    }
  }

  // Get immediate children kd_nodes of the current query and references nodes
  std::vector<ANNkd_node*>* query_nodes = new std::vector<ANNkd_node*>();
  std::vector<ANNkd_node*>* ref_nodes = new std::vector<ANNkd_node*>();
  N_q->child_nodes(*query_nodes), N_r->child_nodes(*ref_nodes);

  // Set parents before recursion
  Rcout << "Recursing" << std::endl;
  N_q_par = N_q, N_r_par = N_r;

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

// The minimum distance between any two points descendent of N_i and N_j
ANNdist DualTree::min_dist(ANNkd_node* N_i, ANNkd_node* N_j){

  Rcout << "Computing the minimum distance" << std::endl;
  std::vector<int>* i_desc_ids = new std::vector<int>();
  std::vector<int>* j_desc_ids = new std::vector<int>();
  N_i->desc_ids(*i_desc_ids);
  N_j->desc_ids(*j_desc_ids);

  ANNdist min_dist = std::numeric_limits<ANNdist>::infinity(), dist;
  for (std::vector<int>::iterator i = i_desc_ids->begin(); i != i_desc_ids->end(); ++i){
    for (std::vector<int>::iterator j = j_desc_ids->begin(); j != j_desc_ids->end(); ++j){
      ANNpoint p_i = this->rtree->pts[*i];
      ANNpoint p_j = this->qtree->pts[*j];
      if ((dist = annDist(d, p_i, p_j)) < min_dist){
        min_dist = dist;
      }
    }
  }
  return min_dist;
}

// The maximum distance between any two points descendent of N_i and N_j
ANNdist DualTree::max_dist(ANNkd_node* N_i, ANNkd_node* N_j){
  std::vector<int>* i_desc_ids = new std::vector<int>();
  std::vector<int>* j_desc_ids = new std::vector<int>();
  N_i->desc_ids(*i_desc_ids);
  N_j->desc_ids(*j_desc_ids);

  ANNdist max_dist = 0, dist;
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

// Get the maximum distance between the centroid of the convex subset (the box of the current node)
// and the points within the children of the current node
ANNdist DualTree::max_child_dist(ANNkd_node* N_i, bool ref_tree){


  std::vector<int> point_ids = std::vector<int>();
  N_i->node_ids(point_ids);

  // Return upper bounds if splitting node
  if (point_ids.size() == 0){
    return this->max_desc_dist(N_i, ref_tree);
  } else {
    ANNpoint centroid = this->centroid(N_i, ref_tree);
    ANNdist max_dist = 0;
    for (std::vector<int>::iterator pt_idx = point_ids.begin(); pt_idx != point_ids.end(); ++pt_idx){
      ANNpoint cpt = ref_tree ? rtree->pts[*pt_idx] : qtree->pts[*pt_idx];
      ANNdist dist = annDist(d, cpt, centroid);
      if (dist > max_dist){ max_dist = dist; }
    }
    return max_dist;
  }
}

// Get the maximum distance between the centroid of the convex subset (the box of the current node)
// and the points within a given node
ANNdist DualTree::max_desc_dist(ANNkd_node* N_i, bool ref_tree){
  ANNpoint centroid = this->centroid(N_i, ref_tree);
  std::vector<int> desc_ids = std::vector<int>();
  N_i->desc_ids(desc_ids);
  ANNdist max_dist = 0;
  for (std::vector<int>::iterator pt_idx = desc_ids.begin(); pt_idx != desc_ids.end(); ++pt_idx){
    ANNpoint cpt = ref_tree ? rtree->pts[*pt_idx] : qtree->pts[*pt_idx];
    ANNdist dist = annDist(d, cpt, centroid);
    if (dist > max_dist){ max_dist = dist; }
  }
  return max_dist;
}


#include "kd_util.h"
#include "ANN/ANNx.h"
#include "ANN/ANNperf.h"

ANNorthRect DualTree::convex_subset(ANNkd_node* N_i, bool ref_tree){
  Rcout << "Getting Convex Subset" << std::endl;
  R_FlushConsole();

  // Choose the tree
  ANNkd_tree* ctree = ref_tree ? rtree : qtree;

  // Get child ids
  std::vector<int> desc_ids = std::vector<int>();
  N_i->desc_ids(desc_ids);
  ANNidxArray pidx = (ANNidxArray) &desc_ids[0];

  // for (int i = 0; i < desc_ids.size(); ++i){
  //   Rcout << pidx[i] << std::endl;
  // }

  // Get the points associated with the child ids
  std::vector<ANNpoint> desc_pts = std::vector<ANNpoint>(desc_ids.size());
  for (int i = 0; i < desc_ids.size(); ++i){ desc_pts.at(i) = (ANNpoint)(ctree->pts[pidx[i]]); }
  ANNpointArray pa = (ANNpointArray) &desc_pts[0];

  // for (int i = 0; i < desc_pts.size(); ++i){
  //   for (int d_i = 0; d_i < d; ++d_i) { Rcout << pa[i][d_i] << ", "; }
  //   Rcout << std::endl;
  // }

  // Compute the bounding box
  ANNorthRect bnd_box(d);
  annEnclRect((ANNpointArray) pa,
              (ANNidxArray) pidx,
              (int) desc_ids.size(),
              (int) d,
              bnd_box); // construct bounding rectangle

  // delete [] pidx;
  // delete [] pa;
  // Return results
  // NumericVector res = NumericVector::create(*bnd_box.lo, *bnd_box.hi);
  // ANNorthRect bnd_box(d);
  return bnd_box;
}

ANNpoint DualTree::centroid(ANNkd_node* N_i, bool ref_tree){
  ANNorthRect bnds = convex_subset(N_i, ref_tree);
  ANNpoint centroid = new ANNcoord[d];
  for (int i = 0; i < d; ++i){
    Rcout << "Bounds (lo, hi): " << bnds.lo[i] << ", " << bnds.hi[i] << std::endl;
    centroid[i] = ANNcoord(bnds.lo[i] + bnds.hi[i] / 2.0);
    Rcout << "Centroid: " << centroid[i] << std::endl;
  }
  return(centroid);
}
