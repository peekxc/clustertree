#include <Rcpp.h>
using namespace Rcpp;

#include "dual_tree.h"
#include "kd_tree.h"

#define NDEBUG // <-- for 'debug' mode
//#undef NDEBUG // <-- for 'production' mode
#ifdef NDEBUG
#include <cassert>
#endif

static int left_count = 0, right_count = 0;
// // static variables to update in the recursion
// #include "pr_queue_k.h"
// ANNmin_k* k_min_pts; // set of k closest points

// DualTree enables the dual-traversal of two kd-tree's at the same time, and
// for computing the associated bounds that come with them
DualTree::DualTree(ANNkd_tree* ref_tree, ANNkd_tree* query_tree) : d(ref_tree->theDim()), knn_identity(ref_tree == query_tree) {
  if (ref_tree->theDim() != query_tree->theDim()){
    stop("Dimensionality of the query set does not match the reference set.");
  }
  //Rcout << "RTree: " <<  ref_tree << "QTree: " << query_tree << std::endl;
  rtree = ref_tree, qtree = query_tree;
  // rtree->Print((ANNbool) true, std::cout);
  //Rcout << "RTree: " <<  rtree << "QTree: " << qtree << std::endl;
  N_q_par = N_r_par = NULL;
  BC_check = new std::map< std::pair<int,int>, bool>();
}

// If the DT's are for KNN, first initialize the appropriate data structures for knn
// This includes a map from a given query pt's index (ANNidx) to a priority queue containing the
// kNN distances (ANNdist) and indices of the corresponing reference points (ANNidx)
void DualTree::knn_initialize(const int k){
  knn = new std::unordered_map<ANNidx, ANNmin_k*>();
  for (int i = 0; i < qtree->n_pts; ++i) {
    knn->insert(std::pair<ANNidx, ANNmin_k*>(qtree->pidx[i], new ANNmin_k(k)));
  }

  // Create map of bounds to memoize the results of the recursion
  bounds = new std::unordered_map<ANNkd_node*, Bound>();
}

#include "kd_util.h"
void DualTree::test_cases(List& res){

  // Print the reference tree
  // rtree->Print((ANNbool) true, std::cout);

  // Initialize KNN structures + bounds map
  knn_initialize(5);

  // Test all methods work for reference tree (recursively)
  Rcout << " ======== Beginning test cases for REFERENCE tree ======== " << std::endl;
  test_cases(res, rtree->root, 0, true);

  // Test all methods work for query tree (recursively)
  // Rcout << " ======== Beginning test cases for QUERY tree ======== " << std::endl;
  // test_cases(res, qtree->root, 0, false);
}

int traversed = 0; // How many recursions have there been?
void DualTree::test_cases(List& res, ANNkd_node* N, int depth, bool ref_tree) {
  traversed++;

  // Get child nodes
  std::vector<ANNkd_node*> child_nodes = std::vector<ANNkd_node*>();
  N->child_nodes(child_nodes);

  // Print the point address to the current node for debugging
  Rcout << "Node ptr: " << N << std::endl;

  // Create bounds object for testing
  if (bounds->find(N) == bounds->end()){ bounds->insert(std::pair<ANNkd_node*, Bound>(N, Bound())); }
  Bound& N_bound = (*bounds)[N];

  // Test bounding box calculation
  ANNorthRect& bnd_box = convex_subset(N, N_bound, ref_tree);
  if (N_bound.bnd_box != NULL){
    Rcout << "Bounding Box: (";
    for (int i = 0; i < d; ++i) { Rcout << N_bound.bnd_box->lo[i] << ", "; }
    Rcout << "), (";
    for (int i = 0; i < d; ++i) { Rcout << N_bound.bnd_box->hi[i] << ", "; }
    Rcout << ")" << std::endl;
  }

  // Test retrieving the centroid
  ANNpoint node_centroid = centroid(N, N_bound, ref_tree);
  if (N_bound.centroid != NULL){
    Rcout << "Centroid: ";
    for (int i = 0; i < d; ++i) { Rcout << N_bound.centroid[i] << ", "; }
    Rcout << std::endl;
  }

  // Test getting the max child distance
  ANNdist n_max_cd = max_child_dist(N, N_bound, ref_tree);
  Rcout << "Max child dist: " << n_max_cd << " (now " << N_bound.rho << ")" << std::endl;

  // Test getting the max desc. distance
  ANNdist n_max_dd = max_desc_dist(N, N_bound, ref_tree);
  Rcout << "Max child dist: " << n_max_dd << " (now " << N_bound.lambda << ")" << std::endl;

  // Recurse
  // Rcout << "Children: ";
  for (std::vector<ANNkd_node*>::iterator node = child_nodes.begin(); node != child_nodes.end(); ++node){
    //Rcout << *node << ", ";
    test_cases(res, (ANNkd_node*) *node, depth + 1, ref_tree);
  }
  // Rcout << std::endl;
//
//   // Test cases: starting at leaves, and working up the tree, test both static and recursive methods
//   Rcout << "Depth: " << depth << ", Traversed: " << traversed << std::endl;
//   R_FlushConsole();
//
//   // Test node id retrieval
//   std::vector<int> node_ids = std::vector<int>();
//   N->node_ids(node_ids);
//   Rcout << "Node ids: ";
//   for (int i = 0; i < node_ids.size(); ++i){ Rcout << node_ids.at(i) << ", "; }
//   Rcout << std::endl;

  // Test Descendent (recursive) node retrieval
  // std::vector<int> desc_ids = std::vector<int>();
  // N->desc_ids(desc_ids);
  // Rcout << "Desc. node ids: ";
  // for (int i = 0; i < desc_ids.size(); ++i){
  //   Rcout << desc_ids.at(i) << ", ";
  // }
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


  // Test bounds checking
  // res["bounds"] = B(N);

  // Return
  return;
}

// Recursive bound on a given query node
ANNdist DualTree::B(ANNkd_node* N_q){
  // If doesn't exist, create a bound object
  if (bounds->find(N_q) == bounds->end()){
    // Rcout << "Inserting Bound object for: " << N_q << std::endl;
    bounds->insert(std::pair<ANNkd_node*, Bound>(N_q, Bound()));
  }

  // If B has been computed before, return it, otherwise look at what needs computed
  Bound& nq_bound = (*bounds)[N_q];
  if (nq_bound.B != -1.0 && nq_bound.B != ANN_DIST_INF){ return(nq_bound.B); }
  else {
    // These bounds need to be computed anyways, so compute them from the start
    max_child_dist(N_q, nq_bound, false);
    max_desc_dist(N_q, nq_bound, false);

    // Rcout << "Computing the bound of a query node (B = " << nq_bound.B << ") " << std::endl;
    ANNdist b1 = 0, b2 = ANN_DIST_INF, b3 = ANN_DIST_INF, b4 = ANN_DIST_INF;

    // Loop through child point ids for bounds 1 and 2
    // std::vector<int> query_ids = std::vector<int>();
    // N_q->node_ids(query_ids);

    ANNdist max_k = bounds->at(N_q).min_knn;
    // Rcout << "Minimax KNN for " << N_q << " : " << max_k << std::endl;
    ANNdist child_bound = max_k + nq_bound.rho + nq_bound.lambda;
    b1 = max_k > b1 ? max_k : b1;
    b2 = child_bound < b2 ? child_bound : b2;

    // for (std::vector<int>::iterator idx = query_ids.begin(); idx != query_ids.end(); ++idx){
    //   assert(nq_bound.rho != 1 && nq_bound.lambda != -1); // rho and lambda should exist at this point
    //   if (knn->find(*idx) != knn->end()){
    //     ANNdist max_k = knn->at(*idx)->max_key();
    //     ANNdist child_bound = max_k + nq_bound.rho + nq_bound.lambda;
    //     b1 = max_k > b1 ? max_k : b1;
    //     b2 = child_bound < b2 ? child_bound : b2;
    //   }
    // }

    // Loop through child kd nodes for bounds 1 and 3
    std::vector<ANNkd_node*> child_nodes = std::vector<ANNkd_node*>();
    N_q->child_nodes(child_nodes);

    if (child_nodes.size() > 0){
      for (std::vector<ANNkd_node*>::iterator N_c = child_nodes.begin(); N_c != child_nodes.end(); ++N_c){
        ANNdist nc_bound = B(*N_c); // Compute child node bound
        assert(bounds->find(*N_c) != bounds->end()); // Child bound should now exist
        ANNdist child_bound = nc_bound + 2 * (nq_bound.lambda - max_desc_dist(*N_c, (*bounds)[*N_c], false));
        b1 = nc_bound > b1 ? nc_bound : b1;
        b3 = child_bound < b3 ? child_bound : b3;
      }
    }

    // ----- Fourth bound -----
    assert(bounds->find(N_q_par) == bounds->end()); // parent should exist at this point
    b4 = N_q_par == NULL ? ANN_DIST_INF : (*bounds)[N_q_par].B;

    // Final bound (lower)
    ANNdist final_bound = std::min((ANNdist) std::min(b1, b2), (ANNdist) std::min(b3, b4));
    nq_bound.B = final_bound;

    // Return final bound
    #ifdef NDEBUG
      // Rcout << N_q << ": final Bound == " << final_bound;
      // Rprintf(" {%.2f, %.2f, %.2f, %.2f, max_cd = %.2f, max_dd = %.2f}\n", b1, b2, b3, b4, nq_bound.rho, nq_bound.lambda);
    #endif
    return final_bound;
  }
}




void DualTree::KNN(int k, NumericMatrix& dists, IntegerMatrix& ids){

  // Initialize map of k-nearest neighbor queues
  // TODO: Map back to integer vector
  knn_initialize(k);

  // Start the search!
  // Rcout << "Starting DFS search!" << std::endl;
  DFS(rtree->root, qtree->root);
  // Rcout << " Total Recursion Calls: " << traversed << std::endl;

  // Copy over the distances and ids
  // Rcout << "Copying distance and ids over to R memory" << std::endl;
  for (int i = 0; i < qtree->n_pts; ++i){
    for (int j = 0; j < k; j++) {		// extract the k-th closest points
      dists(i, j) = knn->at(i)->ith_smallest_key(j);
      ids(i, j) = knn->at(i)->ith_smallest_info(j);
    }
  }

  // Convert to (non-squared) euclidean distances
  for (NumericMatrix::iterator it = dists.begin(); it != dists.end(); *it = sqrt(*it), ++it);

  // deallocate closest point set
  // delete k_min_pts;

  return;
}

// KNN Score function
ANNdist DualTree::Score(ANNkd_node* N_q, ANNkd_node* N_r){
  ANNdist min_dist_qr = min_dist(N_q, N_r);
  // Rcout << "Scoring: Q = " << N_q << ", " <<  "R = " << N_r << std::endl;
  ANNdist bound_nq = B(N_q);
  if (min_dist_qr < bound_nq){
    // Rcout << "Min dist. Q <--> R: " << min_dist_qr << std::endl;
    return (min_dist_qr);
  }
  // Rcout << "Pruning! Q <--> R: " << min_dist_qr << "(> " << bound_nq << ")" << std::endl;
  return ANN_DIST_INF;
}

// key field is distance
// info field is integer (idx)
void DualTree::BaseCase(ANNpoint p_q, ANNpoint p_r, const int q_idx, const int r_idx, ANNkd_node* N_q, ANNkd_node* N_r){
  // Rprintf("dist(%d, %d): ", q_idx, r_idx);
  if (!((bool) ((*BC_check)[std::minmax(q_idx, r_idx)]))){
    ANNdist dist = annDist(d, p_q, p_r);
    // Rprintf("%f \n", dist);
    if (dist < knn->at(q_idx)->max_key()){ knn->at(q_idx)->insert(dist, r_idx); }
    // Rcout << "Storing better neighbor for: " << q_idx << "( id=" << r_idx << ", dist=" << dist << ")" << std::endl;

    // If the query == reference sets, set the other nearest neighbor as well!
    if (knn_identity && q_idx != r_idx && dist < knn->at(r_idx)->max_key()){
      // Rcout << "Storing better neighbor for: " << r_idx << "( id=" << q_idx << ", dist=" << dist << ")" << std::endl;
      knn->at(r_idx)->insert(dist, q_idx);
    }

    // If at any comparison the k-nearest neighbor of any given point is less than the
    // tracked minimax knn of the current query node, update it.
    assert(bounds->find(N_q) != bounds->end());
    if (knn->at(q_idx)->max_key() < bounds->at(N_q).min_knn){
      bounds->at(N_q).min_knn = knn->at(q_idx)->max_key();
      // Rcout << "Updating: Par_q(" << q_idx << ") kNN dist: " << bounds->at(N_q).min_knn << std::endl;
    }
    // If identity, reverse and update other nodes as well
    if (knn_identity && N_q == N_r){
      if (knn->at(r_idx)->max_key() < bounds->at(N_q).min_knn){
        bounds->at(N_q).min_knn = knn->at(r_idx)->max_key();
        // Rcout << "Updating: Par_q(" << r_idx << ") kNN dist: " << bounds->at(N_q).min_knn << std::endl;
      }
    }
  }
  #ifdef NDEBUG
    // else {
    //   Rprintf("<Not Computed>\n");
    // }
  #endif
  (*BC_check)[std::minmax(q_idx, r_idx)] = true; // Make it so other equivalent pairs won't be visited
  // return d;
}

void DualTree::DFS(ANNkd_node* N_q, ANNkd_node* N_r){
  traversed++;
  // Get the points held in both the query and reference nodes
  // Rcout << "Traversing Base cases" << std::endl;
  std::vector<int>* query_pts = new std::vector<int>(), *ref_pts = new std::vector<int>();
  N_q->node_ids(*query_pts), N_r->node_ids(*ref_pts);

  #ifdef NDEBUG
    // if (query_pts->size() > 0){
    //   Rcout << "Query: ";
    //   for (std::vector<int>::iterator q_idx = query_pts->begin(); q_idx != query_pts->end(); ++q_idx){
    //     Rcout << *q_idx << ", ";
    //   }
    //   Rcout << std::endl;
    // }
    //
    // if (ref_pts->size() > 0){
    //   Rcout << "Reference: ";
    //   for (std::vector<int>::iterator r_idx = ref_pts->begin(); r_idx != ref_pts->end(); ++r_idx){
    //     Rcout << *r_idx << ", ";
    //   }
    //   Rcout << std::endl;
    // }
  #endif

  if (Score(N_q, N_r) == ANN_DIST_INF){
    return;
  }

  for (std::vector<int>::iterator q_idx = query_pts->begin(); q_idx != query_pts->end(); ++q_idx){
    for (std::vector<int>::iterator r_idx = ref_pts->begin(); r_idx != ref_pts->end(); ++r_idx){
      BaseCase(qtree->pts[*q_idx], rtree->pts[*r_idx], *q_idx, *r_idx, N_q, N_r); // Pass nodes as well to keep track of min_knn
    }
  }

  // Get immediate children kd_nodes of the current query and references nodes
  std::vector<ANNkd_node*>* query_nodes = new std::vector<ANNkd_node*>();
  std::vector<ANNkd_node*>* ref_nodes = new std::vector<ANNkd_node*>();
  N_q->child_nodes(*query_nodes), N_r->child_nodes(*ref_nodes);

  // Set parents before recursion
  // Rcout << "Setting parents: " << N_q << std::endl;

  // There are at minimum 0, and at maximum 2 nodes per node vector
  int n_query = query_nodes->size(), n_ref = ref_nodes->size();
  if (n_query == 0 && n_ref == 0){ return; }
  else if(n_query == 0 && n_ref > 0) {
    // Traverse through all reference children with the current query node
    N_r_par = N_r; // set current reference node as parent, keep same query parent
    for (std::vector<ANNkd_node*>::iterator n_rc = ref_nodes->begin(); n_rc != ref_nodes->end(); ++n_rc){
      DFS((ANNkd_node*) N_q, (ANNkd_node*) *n_rc);
    }
  } else if (n_query > 0 && n_ref == 0){
    // Traverse through all children
    N_q_par = N_q; // set current query node as parent, keep same reference parent
    for (std::vector<ANNkd_node*>::iterator n_qc = query_nodes->begin(); n_qc != query_nodes->end(); ++n_qc){
      DFS((ANNkd_node*) *n_qc, (ANNkd_node*) N_r);
    }
  } else {
    // Traverse through all children
    N_q_par = N_q, N_r_par = N_r; // set parents as the current nodes
    for (std::vector<ANNkd_node*>::iterator n_qc = query_nodes->begin(); n_qc != query_nodes->end(); ++n_qc){
      for (std::vector<ANNkd_node*>::iterator n_rc = ref_nodes->begin(); n_rc != ref_nodes->end(); ++n_rc){
        DFS((ANNkd_node*) *n_qc, (ANNkd_node*) *n_rc);
      }
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
  if (N_i == N_j) return (ANNdist) 0.0; // Trivial case--nodes are the same
  //Rcout << "Computing the minimum distance" << std::endl;
  std::vector<int>* i_desc_ids = new std::vector<int>();
  std::vector<int>* j_desc_ids = new std::vector<int>();
  N_i->desc_ids(*i_desc_ids);
  N_j->desc_ids(*j_desc_ids);

  ANNdist min_dist = ANN_DIST_INF, dist;
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
      if ((dist = annDist(d, p_i, p_j)) > max_dist){
        max_dist = dist;
      }
    }
  }
  return max_dist;
}

// Get the maximum distance between the centroid of the convex subset (the box of the current node)
// and the points within the children of the current node
ANNdist DualTree::max_child_dist(ANNkd_node* N_i, Bound& ni_bnd, bool ref_tree){

  // If the maximum child distance has already been computed, return it
  if (ni_bnd.rho != -1.0) return ni_bnd.rho;

  // Else look at the points its holding
  std::vector<int> point_ids = std::vector<int>();
  N_i->node_ids(point_ids);

  // Centroid computes (and stores) the bounding box, so access should be fine past here
  ANNpoint ni_centroid = (ni_bnd.centroid == NULL) ? (ANNpoint) this->centroid(N_i, ni_bnd, ref_tree) : ni_bnd.centroid;

  // If it's a splitting node it contains no points. There are two options:
  // 1) Search all of the descendents, computing the true maximum child distance
  // 2) Approximate a looser upper bound using the distance between the precomputed centroid
  // and the bounding rectangle orthogonal to the cutting dimension
  if (point_ids.size() == 0){
    ni_bnd.rho = annDist(d, ni_centroid, ni_bnd.bnd_box->hi); // <-- Looser upper bound
    // Rcout << "== Computed child distance: " << ni_bnd.rho << std::endl;
    // ni_bnd.rho = max_desc_dist(N_i, ni_bnd, ref_tree); // <-- Tighter but slower to compute bound
    return ni_bnd.rho;
  } else {
    // If it's a leaf node and it contains no points, find the one farthest from the centroid
    ANNdist max_dist = 0;
    for (std::vector<int>::iterator pt_idx = point_ids.begin(); pt_idx != point_ids.end(); ++pt_idx){
      ANNpoint cpt = ref_tree ? rtree->pts[*pt_idx] : qtree->pts[*pt_idx];
      ANNdist dist = annDist(d, cpt, ni_centroid);
      // Rcout << "Testing child distance: " << dist << std::endl;
      if (dist > max_dist){ max_dist = dist; }
    }
    ni_bnd.rho = max_dist; // store the bound
    return max_dist;
  }
}

// Get the maximum distance between the centroid of the convex subset (the box of the current node)
// and all of the descendent points within a given node
ANNdist DualTree::max_desc_dist(ANNkd_node* N_i, Bound& ni_bnd, bool ref_tree){

  // If the maximum desc. distance has already been computed, return it
  if (ni_bnd.lambda != -1.0) return ni_bnd.lambda;

  // Get the centroid
  ANNpoint ni_centroid = (ni_bnd.centroid == NULL) ? (ANNpoint) this->centroid(N_i, ni_bnd, ref_tree) : ni_bnd.centroid;

  // Get all descendent point indices
  std::vector<int> desc_ids = std::vector<int>();
  N_i->desc_ids(desc_ids);
  ANNdist max_dist = 0;
  bool check_point_knn = false;
  for (std::vector<int>::iterator pt_idx = desc_ids.begin(); pt_idx != desc_ids.end(); ++pt_idx){
    ANNpoint cpt = ref_tree ? rtree->pts[*pt_idx] : qtree->pts[*pt_idx];
    ANNdist dist = annDist(d, cpt, ni_centroid);
    if (dist > max_dist){ max_dist = dist; }
    // if (knn->at(*pt_idx)->max_key() != ANN_DIST_INF){
    //   check_point_knn = true;
    //   Rcout << "*** " << N_i << " has finite knn points filled!" << std::endl;
    // }
  }
  ni_bnd.lambda = max_dist; // Store the maximum distance by reference for later
  return max_dist;
}


#include "kd_util.h"
#include "ANN/ANNx.h"
#include "ANN/ANNperf.h"

ANNorthRect& DualTree::convex_subset(ANNkd_node* N_i, Bound& ni_bnd, bool ref_tree){
  if (ni_bnd.bnd_box != NULL) return *ni_bnd.bnd_box;
  // Get child ids
  std::vector<int> desc_ids = std::vector<int>();
  N_i->desc_ids(desc_ids);
  ANNidxArray pidx = (ANNidxArray) &desc_ids[0];

  // Compute the bounding rectangle
  ANNorthRect* bnd_box = new ANNorthRect(d); // bounds are dynamically allocated
  annEnclRect((ANNpointArray) ref_tree ? rtree->pts : qtree->pts,
              (ANNidxArray) pidx,
              (int) desc_ids.size(),
              (int) d,
              *bnd_box);
  ni_bnd.bnd_box = bnd_box; // Store bounding rectangle
  return *bnd_box;
}

ANNpoint DualTree::centroid(ANNkd_node* N_i, Bound& ni_bnd, bool ref_tree){
  if (ni_bnd.centroid != NULL) return(ni_bnd.centroid);
  ANNorthRect& bnds = ni_bnd.bnd_box == NULL ? convex_subset(N_i, ni_bnd, ref_tree) : *ni_bnd.bnd_box;
  ANNpoint n_centroid = new ANNcoord[d];
  for (int i = 0; i < d; ++i){
    // Rcout << "Bounds (lo, hi): " << bnds.lo[i] << ", " << bnds.hi[i] << std::endl;
    n_centroid[i] = ANNcoord((bnds.lo[i] + bnds.hi[i]) / 2.0);
    // Rcout << "Centroid: " << centroid[i] << std::endl;
  }
  ni_bnd.centroid = n_centroid;
  return(n_centroid);
}


// [[Rcpp::export]]
List DT_knn(NumericMatrix x, const int k, const int bkt_size = 30) {

  // Copy data over to ANN point array
  int nrow = x.nrow(), ncol = x.ncol();

  ANNpointArray dataPts = annAllocPts(nrow, ncol);
  for(int i = 0; i < nrow; i++){
    for(int j = 0; j < ncol; j++){
      (dataPts[i])[j] = x(i, j);
    }
  }

  // Create kd tree
  ANNkd_tree* kdTree = new ANNkd_tree(dataPts, nrow, ncol, bkt_size, ANN_KD_SUGGEST);

  // Note: the search also returns the point itself (as the first hit)!
  // So we have to look for k+1 points.
  // ANNdistArray dists = new ANNdist[k+1];
  // ANNidxArray nnIdx = new ANNidx[k+1];

  // Distance matrix of kNN distances
  NumericMatrix dists(x.nrow(), k + 1);

  // Id matrix of knn indices
  IntegerMatrix id(x.nrow(), k + 1);

  // Create dual tree using both trees
  DualTree dt = DualTree(kdTree, kdTree);
  List test_res = List();
  //dt.test_cases(test_res);
  dt.KNN(k + 1, dists, id);
  List res = List::create(_["dist"] = dists, _["id"] = id + 1, _["info"] = test_res);

  return(res);
}

/*** R
test_set <- matrix(c(13, 291,
                     57, 145,
                     115, 232,
                     86, 27,
                     145, 28,
                     262, 203,
                     320, 261,
                     379, 174,
                     261, 71,
                     325, 57), byrow=T, ncol=2)
  test_set[, 2] <- 321 - test_set[, 2]

  plot(test_set, xlim = range(test_set[, 1]) + c(-50, 50), ylim = range(test_set[, 2]) + c(-50, 50),
       asp = 1)
  text(test_set, labels = 1:10, pos = 3)
  clustertree:::DT_knn(test_set, k = 4L, bkt_size = 2L)$id[, -1] == unname(dbscan::kNN(test_set, k = 4L)$id)


  clustertree:::DT_knn(test_set, k = 3L)$dist[, -1] == unname(dbscan::kNN(test_set, k = 3L)$dist)


  size <- 1500
  xyz <- as.matrix(data.frame(x = rnorm(size), y = rnorm(size), z = rnorm(size)))
  microbenchmark::microbenchmark(invisible(clustertree:::DT_knn(xyz, k = 30L, bkt_size = 15L)), times = 30)
  microbenchmark::microbenchmark(invisible(dbscan::kNN(xyz, k = 30L)), times = 30)

  */


