#include "dt_knn.h"

ANNkd_node* N_q_par, *N_r_par;

#ifdef ANN_PERF
unsigned int n_traversals = 0;
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
  #ifdef ANN_PERF
    Rcout << "KNN took: " << n_traversals << " traversals\n. Copying to R memory\n";
  #endif
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

  // Compute the minimum box distance
  ANNdist min_dist = annDist(d, nq_bound.centroid, nr_bound.centroid) - nq_bound.lambda - nr_bound.lambda;

  // If it's negative, then the boxes overlap, and the minimum distance between any two points in
  // either branch is 0.
  return(min_dist < 0 ? 0 : min_dist);
}

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

  // If B has been computed before, return it, otherwise look at what needs computed
  BoundKNN& nq_knn_bnd = (BoundKNN&) (*bnd_knn)[N_q];
  if (nq_knn_bnd.B != ANN_DIST_INF){ return(nq_knn_bnd.B); }
  else { // Check to see if we can get better bound

    // Retrieve regular bound object
    const Bound& nq_bound = (const Bound&) (*bounds)[N_q];

    // // There are 4 bounds to compute, although some apply only to leaf or split nodes
    ANNdist bound;
    if (IS_SPLIT(N_q)){
      ANNkd_split* nq_split = AS_SPLIT(N_q); // get split node
      ANNdist lchild_B = (*bnd_knn)[nq_split->child[ANN_LO]].B, rchild_B = (*bnd_knn)[nq_split->child[ANN_HI]].B;
      bound = std::min(
                std::min(lchild_B + 2 * (nq_bound.lambda - (*bounds)[nq_split->child[ANN_LO]].lambda),
                         rchild_B + 2 * (nq_bound.lambda - (*bounds)[nq_split->child[ANN_HI]].lambda)),
                N_q_par == NULL ? ANN_DIST_INF : (*bnd_knn)[N_q_par].B);
    } else {
      ANNkd_leaf* nq_leaf = AS_LEAF(N_q); // get leaf node
      ANNdist max_knn = nq_knn_bnd.maxKNN(nq_leaf->n_pts);
      bound = std::min(std::min(max_knn, nq_knn_bnd.min_knn + nq_bound.rho + nq_bound.lambda),
                       N_q_par == NULL ? ANN_DIST_INF : (*bnd_knn)[N_q_par].B);
    }

    // Final bound (lower)
    nq_knn_bnd.B = bound;
    R_INFO(N_q << ": final Bound == " << nq_knn_bnd.B << (IS_LEAF(N_q) ? " (leaf)" : "(split)"))
    // R_PRINTF(" {%.2f, %.2f, %.2f, %.2f, max_cd = %.2f, max_dd = %.2f}\n", b1, b2, b3, b4, nq_bound.rho, nq_bound.lambda)

    // Return final bound
    return nq_knn_bnd.B;
  }
}

// Simple comparator to sort by value rather than key
template <typename T, typename V>
inline bool sortByValue (const std::pair<T,V>& lhs, const std::pair<T,V>& rhs) { return lhs.second < rhs.second; }

// Pruning dual tree depth-first search traversal
void DualTreeKNN::pDFS(ANNkd_node* N_q, ANNkd_node* N_r) {
  INC_TRAVERSAL(1)

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
  return;
}

// Regular (non-pruning) dual tree DFS traversal
// This shouldn't actually be used in practice, but it's a good function to have to benchmark
// how well the dual traversal performs with and without pruning enabled
void DualTreeKNN::DFS(ANNkd_node* N_q, ANNkd_node* N_r){  INC_TRAVERSAL(1)
  // KD Trees only store points in the leaves; Base case only needed if comparing two leaves
  if (IS_LEAF(N_q) && IS_LEAF(N_r)){ ANN_LEAF(2) // Update leaf counter
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
    ANN_PTS(N_q_leaf->n_pts + N_r_leaf->n_pts)
    return; // If at the base case, don't recurse!
  }

  // --------------------- Begin recursive calls ---------------------
  if (IS_LEAF(N_q) && !IS_LEAF(N_r)){ ANN_LEAF(1); ANN_SPL(1)
    N_r_par = N_r; // set current reference node as parent, keep same query parent
    DFS(N_q, AS_SPLIT(N_r)->child[ANN_LO]); // Go left down the reference node
    DFS(N_q, AS_SPLIT(N_r)->child[ANN_HI]); // Go right down the reference node
  } else if (!IS_LEAF(N_q) && IS_LEAF(N_r)){ ANN_SPL(1); ANN_LEAF(1)
    N_q_par = N_q; // set current query node as parent, keep same reference parent
    DFS(AS_SPLIT(N_q)->child[ANN_LO], N_r); // Go left down the query node
    DFS(AS_SPLIT(N_q)->child[ANN_HI], N_r); // Go right down the query node
  } else { ANN_SPL(4)
    // Set both parents and traverse all combinations
    N_q_par = N_q, N_r_par = N_r;
    DFS(AS_SPLIT(N_q)->child[ANN_LO], AS_SPLIT(N_r)->child[ANN_LO]);
    DFS(AS_SPLIT(N_q)->child[ANN_LO], AS_SPLIT(N_r)->child[ANN_HI]);
    DFS(AS_SPLIT(N_q)->child[ANN_HI], AS_SPLIT(N_r)->child[ANN_LO]);
    DFS(AS_SPLIT(N_q)->child[ANN_HI], AS_SPLIT(N_r)->child[ANN_HI]);
  }
}
