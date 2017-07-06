#include <Rcpp.h>
using namespace Rcpp;

#include "dual_tree.h"
#include "kd_tree.h"

// // static variables to update in the recursion
// #include "pr_queue_k.h"
// ANNmin_k* k_min_pts; // set of k closest points

// DualTree enables the dual-traversal of two kd-tree's at the same time, and
// for computing the associated bounds that come with them
DualTree::DualTree(ANNkd_tree* ref_tree, ANNkd_tree* query_tree) : d(ref_tree->theDim()) {
  if (ref_tree->theDim() != query_tree->theDim()){
    stop("Dimensionality of the query set does not match the reference set.");
  }
  rtree = ref_tree, qtree = query_tree;
  N_q_par = N_r_par = NULL;
}

// If the DT's are for KNN, first initialize the appropriate data structures for knn
// This includes a map from a given query pt's index (ANNidx) to a priority queue containing the
// kNN distances (ANNdist) and indices of the corresponing reference points (ANNidx)
void DualTree::knn_initialize(const int k){
  knn = new std::unordered_map<ANNidx, ANNmin_k*>();
  for (int i = 0; i < qtree->n_pts; ++i) {
    // Rcout << "Inserting: " << qtree->pidx[i] << std::endl;
    knn->insert(std::pair<ANNidx, ANNmin_k*>(qtree->pidx[i], new ANNmin_k(k)));
  }
}

void DualTree::test_cases(List& res){
  // Test all methods work for reference tree (recursively)
  test_cases(res, rtree->root, 0, true);
}
int traversed = 0;
void DualTree::test_cases(List& res, ANNkd_node* N, int depth, bool ref_tree) {
  traversed++;

  // Get child nodes
  std::vector<ANNkd_node*> child_nodes = std::vector<ANNkd_node*>();
  N->child_nodes(child_nodes);

  // Recurse
  for (std::vector<ANNkd_node*>::iterator node = child_nodes.begin(); node != child_nodes.end(); ++node){
    test_cases(res, *node, depth + 1, ref_tree);
  }

  // Test cases: starting at leaves, and working up the tree, test both static and recursive methods
  Rcout << "Depth: " << depth << ", Traversed: " << traversed << std::endl;
  R_FlushConsole();

  // Print the point address to the current node for debugging
  Rcout << "Node ptr: " << N << std::endl;

  // Test node id retrieval
  std::vector<int> node_ids = std::vector<int>();
  N->node_ids(node_ids);
  Rcout << "Node ids: ";
  for (int i = 0; i < node_ids.size(); ++i){
    Rcout << node_ids.at(i) << ", ";
  }
  Rcout << std::endl;

  // Test Descendent (recursive) node retrieval
  std::vector<int> desc_ids = std::vector<int>();
  N->desc_ids(desc_ids);
  Rcout << "Desc. node ids: ";
  for (int i = 0; i < desc_ids.size(); ++i){
    Rcout << desc_ids.at(i) << ", ";
  }
  Rcout << std::endl;

  ANNorthRect bnd_box = convex_subset(N, ref_tree);
  Rcout << "Bounding Box: (";
  for (int i = 0; i < d; ++i) { Rcout << bnd_box.lo[i] << ", "; }
  Rcout << "), (";
  for (int i = 0; i < d; ++i) { Rcout << bnd_box.hi[i] << ", "; }
  Rcout << ")" << std::endl;

  // Test retrieving the centroid
  // ANNpoint node_centroid = centroid(N, ref_tree);
  // Rcout << "Centroid: ";
  // for (int i = 0; i < d; ++i) { Rcout << node_centroid[i] << ", "; }
  // Rcout << std::endl;


  // std::copy(r_centroid, r_centroid + d, r_res.begin());
  // std::copy(q_centroid, q_centroid + d, q_res.begin());
  // res["r_centroid"] = r_res;
  // res["q_centroid"] = q_res;
  //
  // // Test the various node extensions
  // std::vector<int> desc_ids = std::vector<int>(), node_ids = std::vector<int>();
  // std::vector<ANNkd_node*> desc_nodes = std::vector<ANNkd_node*>(), child_nodes = std::vector<ANNkd_node*>();
  // rtree->root->node_ids(node_ids);
  // rtree->root->child_nodes(child_nodes);
  // rtree->root->desc_ids(desc_ids);
  // rtree->root->desc_nodes(desc_nodes);
  //
  // res["node_ids"] = wrap(node_ids); // should be empty
  // res["desc_ids"] = wrap(desc_ids); // should contain all ids, inorder

  // Test KNN structures
  //basicDFS(rtree->root);

  // Test bounds checking
  knn_initialize(5);
  res["bounds"] = B(qtree->root);

  // Return
  return;
}

void DualTree::KNN(int k, NumericMatrix& dists, IntegerMatrix& ids){

  // Initialize map of k-nearest neighbor queues
  // TODO: Map back to integer vector
  knn_initialize(k);

  // Create map of bounds to memoize the results of the recursion
  bounds = new std::unordered_map<ANNkd_node*, Bound>();

  // Start the search!
  Rcout << "Starting DFS search!" << std::endl;
  DFS(rtree->root, qtree->root);

  // Copy over the distances and ids
  Rcout << "Copying distance and ids over to R memory" << std::endl;
  for (int i = 0; i < qtree->n_pts; ++i){
    for (int j = 0; j < k; j++) {		// extract the k-th closest points
      dists(i, j) = knn->at(qtree->pidx[i])->ith_smallest_key(j);
      ids(i, j) = knn->at(qtree->pidx[i])->ith_smallest_info(j);
    }
  }

  // deallocate closest point set
  // delete k_min_pts;

  return;
}

// double DualTree::B1(ANNkd_node* N_q){
//   std::vector<int> child_idx = std::vector<int>();
//   N_q->node_ids(child_idx);
//   double bound = 0, max_k;
//   for (std::vector<int>::iterator idx = child_idx.begin(); idx != child_idx.end(); ++idx){
//     max_k = knn->at(*idx)->max_key();
//     bound = max_k > bound ? max_k : bound;
//   }
//
//   // Get child nodes
//   std::vector<ANNkd_node*> child_nodes = std::vector<ANNkd_node*>();
//   N_q->child_nodes(child_nodes);
//   for (std::vector<ANNkd_node*>::iterator node = child_nodes.begin(); node != child_nodes.end(); ++node){
//     bound = std::max(bound, B1(*node));
//   }
//
//   // Use infinity instead of max double
//   bound = bound == ANN_DIST_INF ? std::numeric_limits<ANNdist>::infinity() : bound;
//
//   return(bound);
// }

// double DualTree::B2(ANNkd_node* N_q){
//   std::vector<int> child_idx = std::vector<int>();
//   N_q->node_ids(child_idx);
//   double bound = std::numeric_limits<ANNdist>::infinity(), min_k;
//   for (std::vector<int>::iterator idx = child_idx.begin(); idx != child_idx.end(); ++idx){
//     min_k = knn->at(*idx)->ith_smallest_key(0);
//     bound = min_k < bound ? min_k : bound;
//   }
//   bound += 2 * max_desc_dist(N_q, false);
//   // Use infinity instead of max double
//   bound = bound == ANN_DIST_INF ? std::numeric_limits<ANNdist>::infinity() : bound;
//
//   return(bound);
// }

// Recursive bound on a given query node
ANNdist DualTree::B(ANNkd_node* N_q){

  // If doesn't exist, create a bound object
  if (bounds->find(N_q) == bounds->end()){
    Rcout << "Inserting Bound object for: " << N_q << std::endl;
    bounds->insert(std::pair<ANNkd_node*, Bound>(N_q, Bound()));
  }

  // If B has been computed before, return it, otherwise look at what needs computed
  Bound nq_bound = bounds->find(N_q)->second;
  if (nq_bound.B != -1.0){
    Rcout << "Key value found: " << bounds->find(N_q)->first << std::endl;
    Rcout << "Using precomputed bound: (" << nq_bound.B << ")" << std::endl;
    return(nq_bound.B);
  } else {
    Rcout << "Computing the bound of a query node" << std::endl;

    ANNdist inf_const = std::numeric_limits<ANNdist>::infinity();
    ANNdist b1 = 0, b2 = inf_const, b3 = inf_const, b4 = inf_const;

    // ----- First bound (upper) -----
    std::vector<int> query_ids = std::vector<int>();
    N_q->node_ids(query_ids);
    for (std::vector<int>::iterator idx = query_ids.begin(); idx != query_ids.end(); ++idx){
      if (knn->find(*idx) != knn->end()){
        ANNdist max_k = knn->at(*idx)->max_key();
        b1 = max_k > b1 ? max_k : b1;
      }
    }

    std::vector<ANNkd_node*> child_nodes = std::vector<ANNkd_node*>();
    N_q->child_nodes(child_nodes);
    if (child_nodes.size() > 0){
      for (std::vector<ANNkd_node*>::iterator idx = child_nodes.begin(); idx != child_nodes.end(); ++idx){
        ANNdist child_bound = B(*idx);
        b1 = child_bound > b1 ? child_bound : b1;
      }
    }
    Rcout << "First bound done" << std::endl;

    // ----- Second bound (lower) -----
    for (std::vector<int>::iterator idx = query_ids.begin(); idx != query_ids.end(); ++idx){
      if (knn->find(*idx) != knn->end()){
        ANNdist child_bound = knn->at(*idx)->max_key();
        ANNdist nq_max_cd = nq_bound.rho == -1.0 ? max_child_dist(N_q, nq_bound, false) : nq_bound.rho;
        ANNdist nq_max_dd = nq_bound.lambda == -1.0 ? max_desc_dist(N_q, nq_bound, false) : nq_bound.lambda;
        child_bound += (nq_max_cd + nq_max_dd);
        b2 = child_bound < b2 ? child_bound : b2;
      }
    }
    Rcout << "Second bound done" << std::endl;

    // ----- Third bound (lower) -----
    for (std::vector<ANNkd_node*>::iterator N_c = child_nodes.begin(); N_c != child_nodes.end(); ++N_c){
      if (bounds->find(*N_c) != bounds->end()){ bounds->insert(std::pair<ANNkd_node*, Bound>(*N_c, Bound())); }
      Bound child_bnd = bounds->find(*N_c)->second;
      ANNdist child_bound = B(*N_c) + 2 * (max_desc_dist(N_q, nq_bound, false) - max_desc_dist(*N_c, child_bnd, false));
      b3 = child_bound < b3 ? child_bound : b3;
    }
    Rcout << "Third bound done" << std::endl;

    // ----- Fourth bound -----
    b4 = N_q_par == NULL ? inf_const : B(N_q_par);
    Rcout << "Fourth bound done" << std::endl;

    // Final bound (lower)
    ANNdist final_bound = std::min(std::min(b1, b2), std::min(b3, b4));
    Rcout << "Final Bound: " << final_bound << std::endl;
    nq_bound.B = final_bound;
    bounds->insert(std::pair<ANNkd_node*, Bound>(N_q, nq_bound));
    return final_bound;
  }
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
  ANNdist min_dist_qr = min_dist(N_q, N_r);
  if (min_dist_qr < B(N_q)){
    return (min_dist_qr);
  }
  return std::numeric_limits<ANNdist>::infinity();
}

// key field is distance
// info field is integer (idx)
ANNdist DualTree::BaseCase(ANNpoint p_q, ANNpoint p_r, const int q_idx, const int r_idx){
  ANNdist dist = annDist(d, p_q, p_r);
  Rcout << "Calling base case for " << q_idx << std::endl;
  if (dist < knn->at(q_idx)->max_key())// TODO: make dynamic programming table for BASECASE pairs
  {
    Rcout << "Storing better neighbor for: " << q_idx << "(" << r_idx << ": " << dist << ")" << std::endl;
    knn->at(q_idx)->insert(dist, r_idx);
  }
  return d;
}

void DualTree::DFS(ANNkd_node* N_q, ANNkd_node* N_r){
  // ANNkd_split* split_query = dynamic_cast<ANNkd_split*>(N_q);
  // ANNkd_split* split_ref = dynamic_cast<ANNkd_split*>(N_r);

  // Get the points held in both the query and reference nodes
  Rcout << "Traversing Base cases" << std::endl;
  std::vector<int>* query_pts = new std::vector<int>();
  std::vector<int>* ref_pts = new std::vector<int>();
  N_q->node_ids(*query_pts), N_r->node_ids(*ref_pts);

  Rcout << "Query: ";
  for (std::vector<int>::iterator q_idx = query_pts->begin(); q_idx != query_pts->end(); ++q_idx){
    Rcout << *q_idx << ", ";
  }
  Rcout << std::endl;

  Rcout << "Reference: ";
  for (std::vector<int>::iterator r_idx = ref_pts->begin(); r_idx != ref_pts->end(); ++r_idx){
    Rcout << *r_idx << ", ";
  }
  Rcout << std::endl;


  if (Score(N_q, N_r) == std::numeric_limits<ANNdist>::infinity()){
    return;
  }

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
  //Rcout << "Recursing" << std::endl;
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

  //Rcout << "Computing the minimum distance" << std::endl;
  std::vector<int>* i_desc_ids = new std::vector<int>();
  std::vector<int>* j_desc_ids = new std::vector<int>();
  N_i->desc_ids(*i_desc_ids);
  N_j->desc_ids(*j_desc_ids);

  ANNdist min_dist = std::numeric_limits<ANNdist>::infinity(), dist;
  for (std::vector<int>::iterator i = i_desc_ids->begin(); i != i_desc_ids->end(); ++i){
    for (std::vector<int>::iterator j = j_desc_ids->begin(); j != j_desc_ids->end(); ++j){
      ANNpoint p_i = this->rtree->pts[*i];
      ANNpoint p_j = this->qtree->pts[*j];
      dist = annDist(d, p_i, p_j);
      if (dist < min_dist){
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
ANNdist DualTree::max_child_dist(ANNkd_node* N_i, Bound& ni_bnd, bool ref_tree){
  if (ni_bnd.rho != 1.0) return ni_bnd.rho;
  std::vector<int> point_ids = std::vector<int>();
  N_i->node_ids(point_ids);
  // Return upper bounds if splitting node
  if (point_ids.size() == 0){
    ni_bnd.rho = max_desc_dist(N_i, ni_bnd, ref_tree);
    return ni_bnd.rho;
  } else {
    ANNpoint ni_centroid = (ni_bnd.centroid == NULL) ? (ANNpoint) this->centroid(N_i, ni_bnd, ref_tree) : ni_bnd.centroid;
    ANNdist max_dist = 0;
    for (std::vector<int>::iterator pt_idx = point_ids.begin(); pt_idx != point_ids.end(); ++pt_idx){
      ANNpoint cpt = ref_tree ? rtree->pts[*pt_idx] : qtree->pts[*pt_idx];
      ANNdist dist = annDist(d, cpt, ni_centroid);
      if (dist > max_dist){ max_dist = dist; }
    }
    ni_bnd.rho = max_dist;
    return max_dist;
  }
}

// Get the maximum distance between the centroid of the convex subset (the box of the current node)
// and the points within a given node
ANNdist DualTree::max_desc_dist(ANNkd_node* N_i, Bound& ni_bnd, bool ref_tree){
  if (ni_bnd.lambda != 1.0) return ni_bnd.lambda;
  ANNpoint ni_centroid = (ni_bnd.centroid == NULL) ? (ANNpoint) this->centroid(N_i, ni_bnd, ref_tree) : ni_bnd.centroid;
  std::vector<int> desc_ids = std::vector<int>();
  N_i->desc_ids(desc_ids);
  ANNdist max_dist = 0;
  for (std::vector<int>::iterator pt_idx = desc_ids.begin(); pt_idx != desc_ids.end(); ++pt_idx){
    ANNpoint cpt = ref_tree ? rtree->pts[*pt_idx] : qtree->pts[*pt_idx];
    ANNdist dist = annDist(d, cpt, ni_centroid);
    if (dist > max_dist){ max_dist = dist; }
  }
  ni_bnd.lambda = max_dist;
  return max_dist;
}


#include "kd_util.h"
#include "ANN/ANNx.h"
#include "ANN/ANNperf.h"

ANNorthRect DualTree::convex_subset(ANNkd_node* N_i, bool ref_tree){
  // Get child ids
  std::vector<int> desc_ids = std::vector<int>();
  N_i->desc_ids(desc_ids);
  ANNidxArray pidx = (ANNidxArray) &desc_ids[0];

  // Compute the bounding rectangle
  ANNorthRect bnd_box(d); // bounds are dynamically allocated
  annEnclRect((ANNpointArray) ref_tree ? rtree->pts : qtree->pts,
              (ANNidxArray) pidx,
              (int) desc_ids.size(),
              (int) d,
              bnd_box);
  return bnd_box;
}

ANNpoint DualTree::centroid(ANNkd_node* N_i, Bound& ni_bnd, bool ref_tree){
  if (ni_bnd.centroid == NULL) return(ni_bnd.centroid);
  ANNorthRect bnds = convex_subset(N_i, ref_tree);
  ANNpoint centroid = new ANNcoord[d];
  for (int i = 0; i < d; ++i){
    // Rcout << "Bounds (lo, hi): " << bnds.lo[i] << ", " << bnds.hi[i] << std::endl;
    centroid[i] = ANNcoord((bnds.lo[i] + bnds.hi[i]) / 2.0);
    // Rcout << "Centroid: " << centroid[i] << std::endl;
  }
  ni_bnd.centroid = centroid;
  return(centroid);
}
