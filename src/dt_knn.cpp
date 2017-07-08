#include <Rcpp.h>
using namespace Rcpp;

#include "dt_knn.h"

// Useful globals to refer to in the recursion
ANNkd_node* N_q, *N_r;
ANNkd_node* N_q_par, *N_r_par;

#ifdef NDEBUG
  unsigned int n_traversals = 0;
  #define INC_TRAVERSAL(x) n_traversals++;
  #define RESET_TRAVERSAL() n_traversals = 0;
#else
  #define INC_TRAVERSAL(x)
  #define RESET_TRAVERSAL()
#endif

// Constructor creates dual tree, detects if reference and query trees are identical
DualTreeKNN::DualTreeKNN(ANNkd_tree* ref_tree, ANNkd_tree* query_tree)
  : DualTree(ref_tree, query_tree), knn_identity(ref_tree == query_tree) {
  N_q_par = N_r_par = NULL;
}

// TODO
// DualTreeKNN::~DualTreeKNN(){
//
// }

// Entry function to start the KNN process. Stores results by reference in dists and ids
void DualTreeKNN::KNN(int k, NumericMatrix& dists, IntegerMatrix& ids, const bool prune) {

  RESET_TRAVERSAL() // make sure traversal count is 0
  use_pruning = prune; // Set pruning

  // Create a map between point indices and their corresponding empty k-nearest neighbors
  knn = new std::unordered_map<ANNidx, ANNmin_k*>();
  for (int i = 0; i < qtree->n_pts; ++i) {
    knn->insert(std::pair<ANNidx, ANNmin_k*>(qtree->pidx[i], new ANNmin_k(k)));
  }

  // If pruning is used, create map of bounds to memoize the results of the recursion
  if (use_pruning){ bounds = new std::unordered_map<ANNkd_node*, BoundKNN>(); }

  // Base case check: given a pair of indices, check whether the base case has been satisfied
  BC_check = new std::map< std::pair<int,int>, bool>();

  // Contains the indices of all of the points descendent of the current node
  // Should not exceed the number of query points. Updated several times recursively.
  qpts = new std::vector<ANNidx>();
  qpts->reserve(qtree->n_pts + 1);

  // If pruning is enabled, use that one. Otherwise use regular DFS w/o computing extra bounds.
  if (prune){ pDFS(rtree->root, qtree->root); } else { DFS(rtree->root, qtree->root); }

  // Copy over the distances and ids
  R_INFO("KNN took: " << n_traversals << " traversals\n. Copying to R memory\n")
  for (int i = 0; i < qtree->n_pts; ++i){
    for (int j = 0; j < k; j++) {		// extract the k-th closest points
      dists(i, j) = knn->at(i)->ith_smallest_key(j);
      ids(i, j) = knn->at(i)->ith_smallest_info(j);
    }
  }
  // Convert to (non-squared) euclidean distances
  transform(dists.begin(), dists.end(), dists.begin(), ::sqrt);

  // Cleanup deallocate closest point set
  // delete k_min_pts;

  return;
}

ANNdist DualTreeKNN::min_dist(ANNkd_node* N_q, ANNkd_node* N_r){ // assume query and reference
  if (N_q == N_r) return 0;
  BoundKNN nq_bound = (*bounds)[N_q], nr_bound = (*bounds)[N_r];
  ANNpoint mu_q = nq_bound.centroid == NULL ? centroid(N_q, nq_bound, false) : nq_bound.centroid;
  ANNpoint mu_r = nr_bound.centroid == NULL ? centroid(N_r, nr_bound, true) : nr_bound.centroid;

  // The lower bound between two nodes requires knowledge of the smallest sphere that will fit in the
  // box contained by a given node
  ANNdist lambda_i = ANN_DIST_INF, lambda_j = ANN_DIST_INF;
  for (int i = 0; i < d; ++i){
    lambda_i = std::min(lambda_i, nq_bound.bnd_box->hi[i] - nq_bound.bnd_box->lo[i]);
    lambda_j = std::min(lambda_j, nr_bound.bnd_box->hi[i] - nr_bound.bnd_box->lo[i]);
  }

  // Return bound
  return(annDist(d, mu_q, mu_r) - lambda_i - lambda_j);
}

// ANNdist DualTreeKNN::max_knn(ANNkd_node* N_q){
//   if (IS_LEAF(N_q)){
//     ANNkd_leaf* leaf_node = dynamic_cast<ANNkd_leaf*>(N_q);
//     ANNidxArray held_pts = leaf_node->bkt;
//     ANNdist max_knn = 0;
//     for (int i = 0; i < leaf_node->n_pts; ++i){
//       if (max_knn < knn->at(leaf_node->bkt[i])->max_key()){
//         max_knn = knn->at(leaf_node->bkt[i])->max_key();
//       }
//     }
//     return max_knn;
//   } else {
//     // ANNkd_leaf
//     // return std::max()
//   }
// };

ANNdist DualTreeKNN::max_knn_B(ANNkd_node* N_q){
  if (IS_LEAF(N_q)) {
    return((*bounds)[N_q].knn_known == AS_LEAF(N_q)->n_pts ? (*bounds)[N_q].max_real_knn : ANN_DIST_INF);
  }
  else {
    ANNkd_split* split_node = AS_SPLIT(N_q);
    return(std::max(max_knn_B(split_node->child[ANN_LO]), max_knn_B(split_node->child[ANN_HI])));
  }
}

// Recursive bound on a given query node
ANNdist DualTreeKNN::B(ANNkd_node* N_q){
  // If doesn't exist, create a bound object
  if (bounds->find(N_q) == bounds->end()){
    R_INFO("Inserting Bound object for: " << N_q << "\n")
    bounds->insert(std::pair<ANNkd_node*, BoundKNN>(N_q, BoundKNN()));
  }

  // If B has been computed before, return it, otherwise look at what needs computed
  BoundKNN& nq_bound = (*bounds)[N_q];
  if (nq_bound.B != -1.0 && nq_bound.B != ANN_DIST_INF){ return(nq_bound.B); }
  else {
    // These bounds need to be computed anyways, so compute them from the start
    max_child_dist(N_q, nq_bound, false); // rho (updated in nq_bound)
    max_desc_dist(N_q, nq_bound, false); // lambda (updated in nq_bound)

    // Get the minimax knn distance for the current query node
    ANNdist max_k = bounds->at(N_q).min_knn;
    // R_INFO("Minimax KNN for " << N_q << " : " << max_k << "\n")

    // This is used in combination with rho and lambda
    ANNdist child_bound = max_k + nq_bound.rho + nq_bound.lambda;

    // Create the various bounds
    ANNdist b1 = max_k, b2 = child_bound, b3 = ANN_DIST_INF, b4 = N_q_par == NULL ? ANN_DIST_INF : (*bounds)[N_q_par].B;

    // Loop through child kd nodes for bounds 1 and 3
    if (IS_SPLIT(N_q)){
      std::vector<ANNkd_node*> child_nodes = std::vector<ANNkd_node*>();
      N_q->child_nodes(child_nodes);
      for (VI(ANNkd_node*) N_c = child_nodes.begin(); N_c != child_nodes.end(); ++N_c){
        ANNdist nc_bound = B(*N_c); // Compute child node bound
        assert(bounds->find(*N_c) != bounds->end()); // Child bound should now exist
        ANNdist child_bound = nc_bound + 2 * (nq_bound.lambda - max_desc_dist(*N_c, (*bounds)[*N_c], false));
        b1 = nc_bound > max_k ? nc_bound : max_k;
        b3 = child_bound < b3 ? child_bound : b3;
      }
    }

    // Final bound (lower)
    nq_bound.B = std::min((ANNdist) std::min(b1, b2), (ANNdist) std::min(b3, b4));
    R_INFO(N_q << ": final Bound == " << nq_bound.B << (IS_LEAF(N_q) ? " (leaf)" : "(split)"))
    R_PRINTF(" {%.2f, %.2f, %.2f, %.2f, max_cd = %.2f, max_dd = %.2f}\n", b1, b2, b3, b4, nq_bound.rho, nq_bound.lambda)

    // Return final bound
    return nq_bound.B;
  }
}

// KNN Score function
// ANNdist DualTreeKNN::Score(ANNkd_node* N_q, ANNkd_node* N_r){
//   R_INFO("Scoring: Q = " << N_q << ", " <<  "R = " << N_r << "\n")
//   ANNdist min_dist_qr = min_dist(N_q, N_r), bound_nq = B(N_q);
//   if (min_dist_qr < bound_nq){
//     R_INFO("Min dist. Q <--> R: " << min_dist_qr << "\n")
//     return (min_dist_qr);
//   }
//   R_INFO("Pruning! Q <--> R: " << min_dist_qr << "(> " << bound_nq << ")" <<  "\n")
//   return ANN_DIST_INF;
// }

//annBoxDistance(q, bnd_box_lo, bnd_box_hi, dim)
// void DualTree::BaseCase(ANNpoint p_q, ANNpoint p_r, const int q_idx, const int r_idx, ANNkd_node* N_q, ANNkd_node* N_r){
//   R_PRINTF("dist(%d, %d): ", q_idx, r_idx);
//   assert(bounds->find(N_q) != bounds->end());
//
//   // Has this pair been considered before?
//   if (!((bool) ((*BC_check)[std::minmax(q_idx, r_idx)]))){
//     ANNdist dist = annDist(d, p_q, p_r);// (ANNdist) ANN_SUM(box_dist, ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));
//     R_PRINTF("%f \n", dist);
//
//     // Is the computed distance less than the current k-nearest neighbor?
//     if (dist < knn->at(q_idx)->max_key()){
//       R_INFO("Storing better neighbor for: " << q_idx << "( id=" << r_idx << ", dist=" << dist << ")\n")
//       knn->at(q_idx)->insert(dist, r_idx);
//     }
//
//     // If the query == reference sets, set the other nearest neighbor as well!
//     if (knn_identity && q_idx != r_idx && dist < knn->at(r_idx)->max_key()){
//       R_INFO("Storing better neighbor for: " << r_idx << "( id=" << q_idx << ", dist=" << dist << ")\n")
//       knn->at(r_idx)->insert(dist, q_idx);
//     }
//     (*BC_check)[std::minmax(q_idx, r_idx)] = true; // Equivalent pairs won't be visited
//   }
// }

// Pruning dual tree depth-first search traversal
void DualTreeKNN::pDFS(ANNkd_node* N_q, ANNkd_node* N_r) {
  // ---------- PRUNING STEP ----------
  // If the score for this node combination is INF, prune the branch
  if (Score(N_q, N_r) == ANN_DIST_INF){
    R_INFO("Combination pruned! (" << N_q << ", " << N_r << ")\n")
    return;
  }
  INC_TRAVERSAL(1)

  // KD Trees only store points in the leaves; Base case only needed if comparing two leaves
  if (IS_LEAF(N_q) && IS_LEAF(N_r)){
    // Get the points held in both the query and reference nodes
    std::vector<int> rpts = std::vector<int>();
    N_r->node_ids(rpts);

    N_q->node_ids(*qpts); // add points
    for (VI(int) q_idx = qpts->begin(); q_idx != qpts->end(); ++q_idx){
      for (VI(int) r_idx = rpts.begin(); r_idx != rpts.end(); ++r_idx){
        R_INFO("Calling base case for: q = " << *q_idx << ", r = " << *r_idx << ")\n")

        // Compute Base case, saving knn ids and distances along the way
        BaseCase(qtree->pts[*q_idx], rtree->pts[*r_idx], *q_idx, *r_idx, N_q); // Pass nodes as well to keep track of min_knn

        // Update bounds where necessary to allow for more pruning in the future
        ANNmin_k& query_knn = *(*knn)[*q_idx], &ref_knn = *(*knn)[*r_idx];

        // If the current k-nearest neighbor distance is less than the current minimum knn
        // distance of the query node, update it. Check also for identical nodes.
        ANNdist min_knn = N_q == N_r ? std::min(ref_knn.max_key(), query_knn.max_key()) : query_knn.max_key();
        if (min_knn < (*bounds)[N_q].min_knn){ (*bounds)[N_q].min_knn = min_knn; }
      }
    }
  }

  // If both are leaves, no need to further recurse
  if (IS_LEAF(N_q) && IS_LEAF(N_r)){
    // Done with this leaf, remove points
    qpts->erase(qpts->end() - AS_LEAF(N_q)->n_pts, qpts->end());
    return;
  }

  // Admittedly an odd choice, multimap should work for prioritizing the recursion
  // since apparently pr_queues don't natively support insertion with priority
  std::multimap<ANNdist, NODE_PAIR, std::less<ANNdist> > score_mm = std::multimap<ANNdist, NODE_PAIR, std::less<ANNdist> >();

  // Get immediate children kd_nodes of the current query and references nodes
  std::vector<ANNkd_node*> qn = std::vector<ANNkd_node*>(), rn = std::vector<ANNkd_node*>();
  N_q->child_nodes(qn), N_r->child_nodes(rn);

  // Reference branches still need to be checked
  // ANNdist c_score; // TODO: benchmark inserting all w/ break vs. checking score per insertion
  if (IS_LEAF(N_q) && !IS_LEAF(N_r)){
    N_r_par = N_r; // set current reference node as parent, keep same query parent
    for (VI(ANNkd_node*) n_rc = rn.begin(); n_rc != rn.end(); ++n_rc){
      score_mm.insert(std::pair<ANNdist, NODE_PAIR>(Score(N_q, *n_rc), NODE_PAIR(N_q, *n_rc)));
    }
  } else if (!IS_LEAF(N_q) && IS_LEAF(N_r)){
    N_q_par = N_q; // set current query node as parent, keep same reference parent
    for (VI(ANNkd_node*) n_qc = qn.begin(); n_qc != qn.end(); ++n_qc){
      score_mm.insert(std::pair<ANNdist, NODE_PAIR>(Score(*n_qc, N_r), NODE_PAIR(*n_qc, N_r)));
    }
  } else {
    N_q_par = N_q, N_r_par = N_r; // set parents as the current nodes
    for (VI(ANNkd_node*) n_qc = qn.begin(); n_qc != qn.end(); ++n_qc){
      for (VI(ANNkd_node*) n_rc = rn.begin(); n_rc != rn.end(); ++n_rc){
        score_mm.insert(std::pair<ANNdist, NODE_PAIR>(Score(*n_qc, *n_rc), NODE_PAIR(*n_qc, *n_rc)));
      }
    }
  }

  // Loop through the map values and recurse combinations with scores < infinity
  // Prioritizes 'closer' branches
  for (std::multimap<ANNdist, NODE_PAIR>::iterator it = score_mm.begin(); it != score_mm.end(); ++it){
    if ((*it).first == ANN_DIST_INF){
      R_INFO("----- Pruning recursion -----\n")
      break; // Don't recurse through INF scored branches
    }
    pDFS(it->second.first, it->second.second);
    // Done with this leaf, remove points
    //qpts->erase(qpts->end() - AS_LEAF(N_q)->n_pts, qpts->end());
  }
  return;
}

// Regular (non-pruning) dual tree DFS traversal
void DualTreeKNN::DFS(ANNkd_node* N_q, ANNkd_node* N_r){
  INC_TRAVERSAL(1)

  // KD Trees only store points in the leaves; Base case only needed if comparing two leaves
  if (IS_LEAF(N_q) && IS_LEAF(N_r)){
    // Get the points held in both the query and reference nodes
    std::vector<int> qpts = std::vector<int>(), rpts = std::vector<int>();
    N_q->node_ids(qpts), N_r->node_ids(rpts);
    for (VI(int) q_idx = qpts.begin(); q_idx != qpts.end(); ++q_idx){
      for (VI(int) r_idx = rpts.begin(); r_idx != rpts.end(); ++r_idx){
        BaseCase(qtree->pts[*q_idx], rtree->pts[*r_idx], *q_idx, *r_idx, N_q); // Pass nodes as well to keep track of min_knn
      }
    }
  }

  // If both are leaves, no need to further recurse
  if (IS_LEAF(N_q) && IS_LEAF(N_r)){ return; }

  // --- This should no longer be needed: each tree should only have two children ---
  // Get immediate children kd_nodes of the current query and references nodes
  // std::vector<ANNkd_node*> qn = std::vector<ANNkd_node*>(), rn = std::vector<ANNkd_node*>();
  // N_q->child_nodes(qn), N_r->child_nodes(rn);

  // Reference branches still need to be checked
  if (IS_LEAF(N_q) && !IS_LEAF(N_r)){
    N_r_par = N_r; // set current reference node as parent, keep same query parent

    // Go left and right
    DFS(N_q, AS_SPLIT(N_r)->child[ANN_LO]);
    DFS(N_q, AS_SPLIT(N_r)->child[ANN_HI]);

    // for (VI(ANNkd_node*) n_rc = rn.begin(); n_rc != rn.end(); ++n_rc){
    //   DFS(N_q, *n_rc);
    // }
  } else if (!IS_LEAF(N_q) && IS_LEAF(N_r)){
    N_q_par = N_q; // set current query node as parent, keep same reference parent
    // Go left and right
    DFS(AS_SPLIT(N_q)->child[ANN_LO], N_r);
    DFS(AS_SPLIT(N_q)->child[ANN_LO], N_r);
    // for (VI(ANNkd_node*) n_qc = qn.begin(); n_qc != qn.end(); ++n_qc){
    //   DFS(*n_qc, N_r);
    // }
  } else {
    N_q_par = N_q, N_r_par = N_r; // set parents as the current nodes

    // Traverse all
    DFS(AS_SPLIT(N_q)->child[ANN_LO], AS_SPLIT(N_r)->child[ANN_LO]);
    DFS(AS_SPLIT(N_q)->child[ANN_LO], AS_SPLIT(N_r)->child[ANN_HI]);
    DFS(AS_SPLIT(N_q)->child[ANN_HI], AS_SPLIT(N_r)->child[ANN_LO]);
    DFS(AS_SPLIT(N_q)->child[ANN_HI], AS_SPLIT(N_r)->child[ANN_HI]);

    // for (VI(ANNkd_node*) n_qc = qn.begin(); n_qc != qn.end(); ++n_qc){
    //   for (VI(ANNkd_node*) n_rc = rn.begin(); n_rc != rn.end(); ++n_rc){
    //     DFS(*n_qc, *n_rc);
    //   }
    // }
  }
}
