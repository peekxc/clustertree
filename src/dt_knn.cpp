#include <Rcpp.h>
using namespace Rcpp;

#include "dt_knn.h"

// Useful globals to refer to in the recursion
ANNkd_node* N_q_par, *N_r_par;

#ifdef NDEBUG
  unsigned int n_traversals = 0;
  #define INC_TRAVERSAL(x) n_traversals++;
  #define RESET_TRAVERSAL() n_traversals = 0;
#else
  #define INC_TRAVERSAL(x)
  #define RESET_TRAVERSAL()
#endif

// Instantiate parent dual tree class along with KNN-specific things, including parent
// node pointers and KNN-bound mapping
DualTreeKNN::DualTreeKNN(const bool prune, const int dim) : DualTree(prune, dim) {
  N_q_par = N_r_par = NULL;
  if (prune){
    R_INFO("Allocating KNN bound map\n")
    bnd_knn = new std::unordered_map<ANNkd_node*, BoundKNN& >();
  }
}

// Now that the tree(s) are made, the knn-specific bound obejcts can be instantiated
void DualTreeKNN::setup(ANNkd_tree* kd_treeQ, ANNkd_tree* kd_treeR) {
  // Call parent class setup to assign trees
  (*this).DualTree::setup(kd_treeQ, kd_treeR);

  R_INFO("Creating KNN bound objects\n")
  // Create new knn bounds using the new kd_node pointers
  if (use_pruning){
    for (std::unordered_map<ANNkd_node*, const Bound& >::iterator bnd = bounds->begin(); bnd != bounds->end(); ++bnd){
      bnd_knn->insert(std::pair<ANNkd_node*, BoundKNN& >(bnd->first, (BoundKNN&) *new BoundKNN()));
    }
  }
}


// TODO
// DualTreeKNN::~DualTreeKNN(){
//
// }

// Entry function to start the KNN process. Stores results by reference in dists and ids
void DualTreeKNN::KNN(int k, NumericMatrix& dists, IntegerMatrix& ids) {
  this->k = k; // set k value
  knn_identity = (qtree == rtree); // is the KNN search using identical trees?

  // Create a map between point indices and their corresponding empty k-nearest neighbors
  BEGIN_PROFILE()
  knn = new std::unordered_map<ANNidx, ANNmin_k*>();
  knn->reserve(qtree->n_pts); // reserve enough buckets to hold at least n_pts items
  for (int i = 0; i < qtree->n_pts; ++i) {
    knn->insert(std::pair<ANNidx, ANNmin_k*>(qtree->pidx[i], new ANNmin_k(k)));
  }
  REPORT_TIME("Creating KNN objects")

  // If pruning is enabled, use that one. Otherwise use regular DFS w/o computing extra bounds.
  RESET_TRAVERSAL() // make sure traversal count is 0
  if (use_pruning){ pDFS(rtree->root, qtree->root); } else { DFS(rtree->root, qtree->root); }

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

  // Assert they exist
  assert(bounds->find(N_q) != bounds->end() && bounds->find(N_q) != bounds->end());

  // Retrieve the bounds
  const Bound& nq_bound = (*bounds)[N_q];
  const Bound& nr_bound = (*bounds)[N_r];

  // // Assert the bounding boxes and centroids were computed
  // assert(nq_bound.bnd_box != NULL && nq_bound.centroid != NULL); // should have bnd_box and centroid precomputed
  // assert(nr_bound.bnd_box != NULL && nr_bound.centroid != NULL); // should have bnd_box and centroid precomputed
  // assert(nq_bound.bnd_box->hi != NULL && nq_bound.bnd_box->lo != NULL);
  // assert(nr_bound.bnd_box->hi != NULL && nr_bound.bnd_box->lo != NULL);
  //
  // // The lower bound between two nodes requires knowledge of the smallest sphere that will fit in the
  // // box contained by a given node. So, find the radius of these spheres.
  // ANNdist lambda_i = ANN_DIST_INF, lambda_j = ANN_DIST_INF;
  // for (int i = 0; i < d; ++i){
  //   lambda_i = std::min(lambda_i, (nq_bound.bnd_box->hi[i] - nq_bound.bnd_box->lo[i])) / 2.0;
  //   lambda_j = std::min(lambda_j, (nr_bound.bnd_box->hi[i] - nr_bound.bnd_box->lo[i])) / 2.0;
  //   // R_INFO("Lambda_i: (" << lambda_i << "), Lambda_j: (" << lambda_j << ")\n")
  // }

  // Return bound
  // R_INFO("NQ centroid: ");
  // for (int i = 0; i < d; ++i) { R_PRINTF("%f ", nq_bound.centroid[i]) }
  // R_INFO("\n")
  // R_INFO("NR centroid: ");
  // for (int i = 0; i < d; ++i) { R_PRINTF("%f ", nr_bound.centroid[i]) }
  // R_INFO("\n")
  ANNdist min_dist = annDist(d, nq_bound.centroid, nr_bound.centroid) - nq_bound.lambda - nr_bound.lambda;
  return(min_dist < 0 ? 0 : min_dist);
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

// Get the maximum KNN distance of any descendent (or child) node
ANNdist DualTreeKNN::max_knn_B(ANNkd_node* N_q){
  // The maximum KNN distance of any node will always be infinite if not even k (potentially non-optimal)
  // neighbors have been found yet. The knn_known property checks this. Note this updates the true maximum
  // knn distance for a given leaf node, and thus, the bound.
  if (IS_LEAF(N_q)) {
    ANNdist max_knn = (*bnd_knn)[N_q].maxKNN(AS_LEAF(N_q)->n_pts);
    (*bnd_knn)[N_q].B = max_knn; // Update leaf maximum knn bound!
    return(max_knn);
  }
  // Otherwise, if dealing with a split node, rather than recurse through the descedents re-checking the max
  // key value, return the maximum  *recorded* KNN distance for the children nodes. Child branches
  // will by default return infinite if they haven't been visited (causing the branch to be recursed into),
  // or if not even k (potentially non-optimal) neighbors of any descendent points have been computed.
  else {
    ANNkd_split* split_node = AS_SPLIT(N_q);
    ANNdist max_knn = std::max((*bnd_knn)[split_node->child[ANN_LO]].B,
                               (*bnd_knn)[split_node->child[ANN_HI]].B);
    (*bnd_knn)[N_q].B = max_knn; // Update leaf maximum knn bound!
    return(max_knn);
  }
}

// Recursive bound on a given query node
ANNdist DualTreeKNN::B(ANNkd_node* N_q){
  // If doesn't exist, create a bound object
  // if (bounds->find(N_q) == bounds->end()){
  //   R_INFO("Inserting Bound object for: " << N_q << "\n")
  //   bounds->insert(std::pair<ANNkd_node*, Bound& >(N_q, (Bound&) *new Bound()));
  // }

  // If B has been computed before, return it, otherwise look at what needs computed
  BoundKNN& nq_knn_bnd = (BoundKNN&) (*bnd_knn)[N_q];
  if (nq_knn_bnd.B != -1.0 && nq_knn_bnd.B != ANN_DIST_INF){ return(nq_knn_bnd.B); }
  else {

    // Retrieve regular bound object
    const Bound& nq_bound = (const Bound&) (*bounds)[N_q];

    // There are 4 bounds
    ANNdist b1, b2, b3, b4;
    if (IS_SPLIT(N_q)){
      ANNkd_split* nq_split = AS_SPLIT(N_q); // get split node
      b1 = b2 = ANN_DIST_INF; // Bounds 1 and 2 don't need to be computed if split node; doesn't hold any points!
      ANNdist lchild_B = B(nq_split->child[ANN_LO]), rchild_B = B(nq_split->child[ANN_HI]);
      b3 = std::min(lchild_B + 2 * (nq_bound.lambda - (*bounds)[nq_split->child[ANN_LO]].lambda),
                    rchild_B + 2 * (nq_bound.lambda - (*bounds)[nq_split->child[ANN_HI]].lambda));
      b4 = B(N_q_par);
    } else {
      ANNkd_leaf* nq_leaf = AS_LEAF(N_q); // get leaf node
      ANNdist max_knn = (nq_knn_bnd.knn_known < nq_leaf->n_pts) ? ANN_DIST_INF : nq_knn_bnd.max_real_knn;
      b1 = max_knn;
      b2 = ANN_DIST_INF;
      for (int i = 0; i < nq_leaf->n_pts; ++i){
        ANNdist b2_tmp = knn->at(nq_leaf->bkt[i])->max_key() + nq_bound.rho + nq_bound.lambda;
        b2 = b2_tmp < b2 ? b2_tmp : b2;
      }
      b3 = ANN_DIST_INF; // Leaf has no children
      b4 = B(N_q_par);
    }

    // ANNdist max_k = IS_SPLIT(N_q) ? ANN_DIST_INF : (nq_knn_bnd.knn_known < AS_LEAF(N_q)->n_pts) ? ANN_DIST_INF : nq_knn_bnd.max_real_knn;
    // ANNdist b2 = nq_knn_bnd.max_real_knn + nq_bound.rho + nq_bound.lambda;
    // R_INFO("Minimax KNN for " << N_q << " : " << max_k << "\n")
    //
    // // Create the various bounds
    // ANNdist b1 = max_k, b3 = ANN_DIST_INF, b4 = N_q_par == NULL ? ANN_DIST_INF : (*bnd_knn)[N_q_par].B;
    //
    // if (IS_SPLIT(N_q)){
    //   ANNkd_split* nq_split = AS_SPLIT(N_q);
    //   ANNdist lchild_B = B(nq_split->child[ANN_LO]), rchild_B = B(nq_split->child[ANN_HI]);
    //   nc_bound + 2 * (nq_bound.lambda - (*bounds)[nq_split->child[ANN_LO]].lambda);
    // }

    // Loop through child kd nodes for bounds 1 and 3
    // if (IS_SPLIT(N_q)){
    //   std::vector<ANNkd_node*> child_nodes = std::vector<ANNkd_node*>();
    //   N_q->child_nodes(child_nodes);
    //   for (VI(ANNkd_node*) N_c = child_nodes.begin(); N_c != child_nodes.end(); ++N_c){
    //     ANNdist nc_bound = B(*N_c); // Compute child node bound
    //     assert(bounds->find(*N_c) != bounds->end()); // Child bound should now exist
    //     ANNdist child_bound = nc_bound + 2 * (nq_bound.lambda - max_desc_dist(*N_c, (*bounds)[*N_c], false));
    //     b1 = nc_bound > max_k ? nc_bound : max_k;
    //     b3 = child_bound < b3 ? child_bound : b3;
    //   }
    // }

    // Final bound (lower)
    nq_knn_bnd.B = std::min((ANNdist) std::min(b1, b2), (ANNdist) std::min(b3, b4));
    R_INFO(N_q << ": final Bound == " << nq_knn_bnd.B << (IS_LEAF(N_q) ? " (leaf)" : "(split)"))
    R_PRINTF(" {%.2f, %.2f, %.2f, %.2f, max_cd = %.2f, max_dd = %.2f}\n", b1, b2, b3, b4, nq_bound.rho, nq_bound.lambda)

    // Return final bound
    return nq_knn_bnd.B;
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


// Simple comparator to sort by value rather than key
template <typename T, typename V>
inline bool sortByValue (const std::pair<T,V>& lhs, const std::pair<T,V>& rhs) { return lhs.second < rhs.second; }

// Pruning dual tree depth-first search traversal
void DualTreeKNN::pDFS(ANNkd_node* N_q, ANNkd_node* N_r) {
  INC_TRAVERSAL(1)
// #ifdef ANN_PERF
//   annUpdateStats();
// #endif

  // KD Trees only store points in the leaves; Base case only needed if comparing two leaves
  if (IS_LEAF(N_q) && IS_LEAF(N_r)) {
    ANN_LEAF(2) // Update leaf counter
    BaseCaseIdentity(AS_LEAF(N_q), AS_LEAF(N_r));
    return; // If at the base case, don't recurse!
  }

  // --------------- Begin prioritized recursive calls ---------------
  if (IS_LEAF(N_q) && !IS_LEAF(N_r)){ ANN_LEAF(1); ANN_SPL(1)
    N_r_par = (ANNkd_node*) N_r; // set current reference node as parent, keep same query parent
    const ANNkd_split* nr_spl = AS_SPLIT(N_r);

    // Sort closer branches
    NodeScore scores[2] = { NodeScore(N_q, nr_spl->child[ANN_LO], Score(N_q, nr_spl->child[ANN_LO])),
                            NodeScore(N_q, nr_spl->child[ANN_HI], Score(N_q, nr_spl->child[ANN_HI])) };
    sort2_sn_stable(scores); // Sort by score using sorting network

    // Visit closer branches first, return at first sign of infinity
    R_PRINTF("Score values: (%c, %c) => %f, (%c, %c) => %f\n",
             node_labels.at(scores[0].lhs), node_labels.at(scores[0].rhs), scores[0].score,
             node_labels.at(scores[1].lhs), node_labels.at(scores[1].rhs), scores[1].score)
    if (scores[0].score == ANN_DIST_INF) return; else pDFS(scores[0].lhs, scores[0].rhs);
    if (scores[1].score == ANN_DIST_INF) return; else pDFS(scores[1].lhs, scores[1].rhs);
  } else if (!IS_LEAF(N_q) && IS_LEAF(N_r)){ ANN_SPL(1); ANN_LEAF(1)
    N_q_par = N_q; // set current query node as parent, keep same reference parent
    ANNkd_split *const nq_spl = AS_SPLIT(N_q);

    // Sort closer branches
    NodeScore scores[2] = { NodeScore(nq_spl->child[ANN_LO], N_r, Score(nq_spl->child[ANN_LO], N_r)),
                            NodeScore(nq_spl->child[ANN_HI], N_r, Score(nq_spl->child[ANN_HI], N_r)) };
    sort2_sn_stable(scores); // Sort by score using sorting network

    // Visit closer branches first, return at first sign of infinity
    R_PRINTF("Score values: (%c, %c) => %f, (%c, %c) => %f\n",
             node_labels.at(scores[0].lhs), node_labels.at(scores[0].rhs), scores[0].score,
             node_labels.at(scores[1].lhs), node_labels.at(scores[1].rhs), scores[1].score)
    if (scores[0].score == ANN_DIST_INF) return; else pDFS(scores[0].lhs, scores[0].rhs);
    if (scores[1].score == ANN_DIST_INF) return; else pDFS(scores[1].lhs, scores[1].rhs);
  } else { ANN_SPL(2)
    N_q_par = N_q, N_r_par = N_r; // set parents as the current nodes

    ANNkd_split *const nq_spl = AS_SPLIT(N_q), *const nr_spl = AS_SPLIT(N_r);
    if (N_q != N_r){
      NodeScore scores[4] = { NodeScore(nq_spl->child[ANN_LO], nr_spl->child[ANN_LO], Score(nq_spl->child[ANN_LO], nr_spl->child[ANN_LO])),
                              NodeScore(nq_spl->child[ANN_LO], nr_spl->child[ANN_HI], Score(nq_spl->child[ANN_LO], nr_spl->child[ANN_HI])),
                              NodeScore(nq_spl->child[ANN_HI], nr_spl->child[ANN_LO], Score(nq_spl->child[ANN_HI], nr_spl->child[ANN_LO])),
                              NodeScore(nq_spl->child[ANN_HI], nr_spl->child[ANN_HI], Score(nq_spl->child[ANN_HI], nr_spl->child[ANN_HI])) };
      sort4_sn_stable(scores); // Sort by score using sorting network

      // Visit closer branches first, return at first sign of infinity
      R_PRINTF("Score values: (%c, %c) => %f, (%c, %c) => %f, (%c, %c) => %f, (%c, %c) => %f\n",
               node_labels.at(scores[0].lhs), node_labels.at(scores[0].rhs), scores[0].score,
               node_labels.at(scores[1].lhs), node_labels.at(scores[1].rhs), scores[1].score,
               node_labels.at(scores[2].lhs), node_labels.at(scores[2].rhs), scores[2].score,
               node_labels.at(scores[3].lhs), node_labels.at(scores[3].rhs), scores[3].score)
      if (scores[0].score == ANN_DIST_INF) return; else pDFS(scores[0].lhs, scores[0].rhs);
      if (scores[1].score == ANN_DIST_INF) return; else pDFS(scores[1].lhs, scores[1].rhs);
      if (scores[2].score == ANN_DIST_INF) return; else pDFS(scores[2].lhs, scores[2].rhs);
      if (scores[3].score == ANN_DIST_INF) return; else pDFS(scores[3].lhs, scores[3].rhs);
    } else {
      // If the nodes are identical, it is unnecessary to score and recurse into both children
      NodeScore scores[3] = { NodeScore(nq_spl->child[ANN_LO], nr_spl->child[ANN_LO], Score(nq_spl->child[ANN_LO], nr_spl->child[ANN_LO])),
                              NodeScore(nq_spl->child[ANN_LO], nr_spl->child[ANN_HI], Score(nq_spl->child[ANN_LO], nr_spl->child[ANN_HI])),
                              NodeScore(nq_spl->child[ANN_HI], nr_spl->child[ANN_HI], Score(nq_spl->child[ANN_HI], nr_spl->child[ANN_HI])) };
      sort3_sn_stable(scores); // Sort by score using sorting network

      // Visit closer branches first, return at first sign of infinity
      R_PRINTF("Score values: (%c, %c) => %f, (%c, %c) => %f, (%c, %c) => %f\n",
               node_labels.at(scores[0].lhs), node_labels.at(scores[0].rhs), scores[0].score,
               node_labels.at(scores[1].lhs), node_labels.at(scores[1].rhs), scores[1].score,
               node_labels.at(scores[2].lhs), node_labels.at(scores[2].rhs), scores[2].score)
      if (scores[0].score == ANN_DIST_INF) return; else pDFS(scores[0].lhs, scores[0].rhs);
      if (scores[1].score == ANN_DIST_INF) return; else pDFS(scores[1].lhs, scores[1].rhs);
      if (scores[2].score == ANN_DIST_INF) return; else pDFS(scores[2].lhs, scores[2].rhs);
    }
  }
//
//
//     nodes[0] = NODE_PAIR(N_q_split->child[ANN_LO], N_r_split->child[ANN_LO]);
//     scores[0] = std::make_pair(0, Score(nodes[0].first, nodes[0].second));
//
//     // R_INFO("Computing score 2 (" << N_q << ", " << N_r_split->child[1] << ")\n")
//     nodes[1] = NODE_PAIR(N_q_split->child[ANN_LO], N_r_split->child[ANN_HI]);
//     scores[1] = std::make_pair(1, Score(nodes[1].first, nodes[1].second));
//
//
//     if (!knn_identity || N_q != N_r){
//       //R_INFO("Computing score 3 (" << N_q_split->child[0] << ", " << N_r << ")\n")
//       nodes[2] = NODE_PAIR(N_q_split->child[ANN_HI], N_r_split->child[ANN_LO]);
//       scores[2] = std::make_pair(2, Score(nodes[2].first, nodes[2].second));
//     } else { scores[2] = std::make_pair(2, ANN_DIST_INF); }
//
//     //R_INFO("Computing score 4 (" << N_q_split->child[1] << ", " << N_r << ")\n")
//     nodes[3] = NODE_PAIR(N_q_split->child[ANN_HI], N_r_split->child[ANN_HI]);
//     scores[3] = std::make_pair(3, Score(nodes[3].first, nodes[3].second));
//   }
//
//   // Loop through the map values and recurse combinations with scores < infinity
//   // Prioritizes 'closer' branches, where 'closer' is determined by the minimum distance between
//   // nodes (and their desc. points). These bounds are precomputed during the tree construction.
//   // std::sort(pr.scores.begin(), pr.scores.end());
//   std::sort(std::begin(scores), std::end(scores), sortByValue<ANNidx, ANNdist>); // sorts by value
//   R_PRINTF("Score idx and values: %d %d %d %d : %f %f %f %f\n", scores[0].first, scores[1].first, scores[2].first,
//   scores[3].first, scores[0].second, scores[1].second, scores[2].second, scores[3].second)
//   // R_PRINTF("[0] = (%p, %p) --> (s: %f)\n", nodes[scores[0].first].first, nodes[scores[0].first].second, scores[0].second)
//   // R_PRINTF("[1] = (%p, %p) --> (s: %f)\n", nodes[scores[1].first].first, nodes[scores[1].first].second, scores[1].second)
//   // R_PRINTF("[2] = (%p, %p) --> (s: %f)\n", nodes[scores[2].first].first, nodes[scores[2].first].second, scores[2].second)
//   // R_PRINTF("[3] = (%p, %p) --> (s: %f)\n", nodes[scores[3].first].first, nodes[scores[3].first].second, scores[3].second)
//   for (int i = 0, idx; i < 4; ++i){
//     if (scores[i].second == ANN_DIST_INF){
//       R_INFO("---- Pruning further recursions ----\n")
//       break;
//     }
//     idx = scores[i].first;
//     pDFS(nodes[idx].first, nodes[idx].second);
//   }
  return;
}

// Regular (non-pruning) dual tree DFS traversal
void DualTreeKNN::DFS(ANNkd_node* N_q, ANNkd_node* N_r){
  INC_TRAVERSAL(1)

  // KD Trees only store points in the leaves; Base case only needed if comparing two leaves
  // BEGIN_PROFILE()
  if (IS_LEAF(N_q) && IS_LEAF(N_r)){
    ANN_LEAF(2) // Update leaf counter
    ANNkd_leaf* N_q_leaf = AS_LEAF(N_q), *N_r_leaf = AS_LEAF(N_r);
    for (int q_i = 0; q_i < N_q_leaf->n_pts; ++q_i){
      for (int r_i = 0; r_i < N_r_leaf->n_pts; ++r_i){
        int q_idx = N_q_leaf->bkt[q_i], r_idx = N_r_leaf->bkt[r_i];

        // Compute Base case, saving knn ids and distances along the way
        if (!((bool) ((*BC_check)[std::minmax(q_idx, r_idx)]))){ // Has this pair been considered before?
          R_INFO("Calling base case for: q = " << q_idx << ", r = " << r_idx << ")\n")
          BaseCase(qtree->pts[q_idx], rtree->pts[r_idx], q_idx, r_idx, N_q); // Pass nodes as well to keep track of min_knn
          (*BC_check)[std::minmax(q_idx, r_idx)] = true; // Equivalent pairs won't be visited again
        }
      }
    }
    return; // If at the base case, don't recurse!
  }
  // REPORT_TIME("Running base cases")

  // -----------------------------------------------------------------
  // --------------------- Begin recursive calls ---------------------
  // -----------------------------------------------------------------
  if (IS_LEAF(N_q) && !IS_LEAF(N_r)){
    ANN_LEAF(1); ANN_SPL(1)
    N_r_par = N_r; // set current reference node as parent, keep same query parent
    DFS(N_q, AS_SPLIT(N_r)->child[ANN_LO]); // Go left down the reference node
    DFS(N_q, AS_SPLIT(N_r)->child[ANN_HI]); // Go right down the reference node
  } else if (!IS_LEAF(N_q) && IS_LEAF(N_r)){
    ANN_SPL(1); ANN_LEAF(1)
    N_q_par = N_q; // set current query node as parent, keep same reference parent
    DFS(AS_SPLIT(N_q)->child[ANN_LO], N_r); // Go left down the query node
    DFS(AS_SPLIT(N_q)->child[ANN_HI], N_r); // Go right down the query node
  } else {
    // Set both parents and traverse all combinations
    N_q_par = N_q, N_r_par = N_r;
    ANN_SPL(2)
    DFS(AS_SPLIT(N_q)->child[ANN_LO], AS_SPLIT(N_r)->child[ANN_LO]);
    DFS(AS_SPLIT(N_q)->child[ANN_LO], AS_SPLIT(N_r)->child[ANN_HI]);
    DFS(AS_SPLIT(N_q)->child[ANN_HI], AS_SPLIT(N_r)->child[ANN_LO]);
    DFS(AS_SPLIT(N_q)->child[ANN_HI], AS_SPLIT(N_r)->child[ANN_HI]);
  }
}
