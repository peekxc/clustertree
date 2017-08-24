#include "dt_knn.h"

#ifdef ANN_PERF
unsigned int n_traversals = 0;
#endif

// Instantiate parent dual tree class along with KNN-specific things, including parent
// node pointers and KNN-bound mapping
DualTreeKNN::DualTreeKNN(const NumericMatrix& q_x, Metric& m, NumericMatrix& r_x, List& config) : DualTree(q_x, m, r_x, config) {
  N_q_par = N_r_par = NULL; // Set parent nodes as null to prepare for traversal
  bnd_knn = std::unordered_map<ANNkd_node*, BoundKNN& >();
  bnd_knn.reserve(bounds.size());

  // Create an initial map between the node pointers and their subsequent bounds (will have to anyways)
  for (std::unordered_map<ANNkd_node*, const Bound& >::iterator bnd = bounds.begin(); bnd != bounds.end(); ++bnd){
    BoundKNN& new_bnd = *new BoundKNN();
    bnd_knn.insert(std::pair<ANNkd_node*, BoundKNN&>((*bnd).first, new_bnd));
  }
}

// Destructor
DualTreeKNN::~DualTreeKNN(){
  N_q_par = NULL;
  N_r_par = NULL;
  k = 0;
  std::unordered_map<ANNkd_node*, BoundKNN& > empty_bnds;
  std::swap(bnd_knn, empty_bnds);
}

// Entry function to start the KNN process. Stores results by reference in dists and ids
List DualTreeKNN::KNN(int k) {
  this->k = k; // set k value
  const int n = qtree->n_pts; // number of query points

  // Create a map between point indices and their corresponding empty k-nearest neighbors
  knn = std::unordered_map<ANNidx, ANNmin_k*>();
  knn.reserve(n); // reserve enough buckets to hold at least n_pts items
  // D = new std::vector<ANNdist>(qtree->n_pts, ANN_DIST_INF); // cache of knn distances
  for (int i = 0; i < n; ++i) {
    knn.insert(std::pair<ANNidx, ANNmin_k*>(qtree->pidx[i], new ANNmin_k(k)));
  }

  // If pruning is enabled, use that one. Otherwise use regular DFS w/o computing extra bounds.
  if (use_pruning){ pDFS(rtree->root, qtree->root); } else { DFS(rtree->root, qtree->root); }

  NumericMatrix dists(n, k - 1);   // Distance matrix of kNN distances
  IntegerMatrix id(n, k - 1);  // Id matrix of knn indices
  NumericVector dist_row = NumericVector(k); // The current distance row
  IntegerVector id_row = IntegerVector(k); // the id row
  for (int i = 0; i < qtree->n_pts; ++i){
    for (int j = 0; j < k; j++) {		// extract the k-th closest points
      dist_row(j) = knn.at(i)->ith_smallest_key(j);
      id_row(j) = knn.at(i)->ith_smallest_info(j);
    }
    LogicalVector take = id_row != i;
    dists(i, _) = as<NumericVector>(dist_row[take]);
    id(i, _) = as<IntegerVector>(id_row[take]);
  }

  // Finalize results for R-end
  m_dist.finalize(dists); // apply final transformation
  id = id + 1; // switch to 1-based for R

  // Cleanup KNN stuff
  std::unordered_map<ANNidx, ANNmin_k*> empty_knn;
  std::swap(knn, empty_knn);

  // Return
  List res = List::create(_["dist"] = dists, _["id"] = id);
  return res;
}

// Reset the knn pririty queues
void DualTreeKNN::resetKNN(){
  std::unordered_map<ANNidx, ANNmin_k*>().swap(knn); // clear old knn
  knn = std::unordered_map<ANNidx, ANNmin_k*>();
  knn.reserve(qtree->n_pts); // reserve enough buckets to hold at least n_pts items
  for (int i = 0; i < qtree->n_pts; ++i) {
    knn.insert(std::pair<ANNidx, ANNmin_k*>(qtree->pidx[i], new ANNmin_k(k)));
  }
}

// Get the maximum KNN distance of any descendent (or child) node
ANNdist DualTreeKNN::max_knn_B(ANNkd_node* N_q){
  // The maximum KNN distance of any node will always be infinite if not even k (potentially non-optimal)
  // neighbors have been found yet. The knn_known property checks this. Note this updates the true maximum
  // knn distance for a given leaf node, and thus, the bound.
  if (IS_LEAF(N_q)) {
    ANNdist max_knn = bnd_knn[N_q].maxKNN(AS_LEAF(N_q)->n_pts);
    bnd_knn[N_q].B = max_knn; // Update leaf maximum knn bound!
    return(max_knn);
  }
  // Otherwise, if dealing with a split node, rather than recurse through the descedents re-checking the max
  // key value, return the maximum  *recorded* KNN distance for the children nodes. Child branches
  // will by default return infinite if they haven't been visited (causing the branch to be recursed into),
  // or if not even k (potentially non-optimal) neighbors of any descendent points have been computed.
  else {
    ANNkd_split* split_node = AS_SPLIT(N_q);
    ANNdist max_knn = std::max(bnd_knn[split_node->child[ANN_LO]].B,
                               bnd_knn[split_node->child[ANN_HI]].B);
    bnd_knn[N_q].B = max_knn; // Update leaf maximum knn bound!
    return(max_knn);
  }
}

// Reset the KNN-related bounds
void DualTreeKNN::resetBounds(){
  for (std::unordered_map<ANNkd_node*, BoundKNN& >::iterator it = bnd_knn.begin(); it != bnd_knn.end(); ++it){
    bnd_knn.at((*it).first).reset();
  }
}

// Update the various properties of the KNN-only bounds related to the DT framework.
// This includes:
// 1. The node-wide number of non-inf (known) KNN distances
// 2. The node-wide minimum kth-nearest neighbor distance
// 3. The node-wide maximum non-inf kth-nearest neighbor distance
// TODO: FIX THE INFORMATION OUTPut AND GET BACK TO FIXED KNN AND DTB
ANNdist DualTreeKNN::updateBounds(ANNdist new_dist, ANNkd_node* N_q_leaf, ANNkd_node* N_r_leaf,
                                  ANNdist min_dist_q, ANNdist min_dist_r,
                                  bool add_known_query, bool add_known_ref){
  // Retrieve Bound objects
  BoundKNN& qbnd = bnd_knn.at(N_q_leaf), &rbnd = bnd_knn.at(N_r_leaf);

  // 'new_dist' encodes a newly calculated distance between q and r. If this distance is non-inf,
  // and the old currently cached distance is inf, then update the number of points known by the
  // current node to have non-inf distances
  if (add_known_query) qbnd.knn_known++;
  if (add_known_ref) rbnd.knn_known++;

  // Update minimum kth nearest neighbor distances
  if (min_dist_q < qbnd.min_knn) qbnd.min_knn = min_dist_q; // Update minimum kth nearest distance on query node
  if (min_dist_r < rbnd.min_knn) rbnd.min_knn = min_dist_r; // Update minimum kth nearest distance on reference node

  // Update maximum (non-inf) distances as well
  if (min_dist_q > qbnd.max_real_knn && min_dist_q != ANN_DIST_INF) qbnd.max_real_knn = min_dist_q;
  if (min_dist_r > rbnd.max_real_knn && min_dist_r != ANN_DIST_INF) rbnd.max_real_knn = min_dist_r;

  // Informational output
  R_INFO("Q (" << node_labels.at(N_q_leaf) << ") max (non-inf) knn: " << qbnd.max_real_knn
               << " (" << qbnd.knn_known << "/" << AS_LEAF(N_q_leaf)->n_pts << " known)"
               << " max_knn: " << qbnd.maxKNN(AS_LEAF(N_q_leaf)->n_pts)
               << " min_knn: " << qbnd.min_knn << "\n")
  R_INFO("R (" << node_labels.at(N_r_leaf) << ") max (non-inf) knn: " << rbnd.max_real_knn
               << " (" << rbnd.knn_known << "/" << AS_LEAF(N_r_leaf)->n_pts << " known)"
               << " max_knn: " << rbnd.maxKNN(AS_LEAF(N_r_leaf)->n_pts)
               << " min_knn: " << rbnd.min_knn << "\n")
  return(new_dist);
}


// Recursive bound on a given query node
ANNdist DualTreeKNN::B(ANNkd_node* N_q){

  // If B has been computed before, return it, otherwise look at what needs computed
  BoundKNN& nq_knn_bnd = bnd_knn[N_q];
  if (nq_knn_bnd.B != ANN_DIST_INF){ return(nq_knn_bnd.B); }
  else { // Check to see if we can get better bound

    // Retrieve regular bound object
    const Bound& nq_bound = bounds[N_q];

    // // There are 4 bounds to compute, although some apply only to leaf or split nodes
    ANNdist bound;
    if (IS_SPLIT(N_q)){
      ANNkd_split* nq_split = AS_SPLIT(N_q); // get split node
      ANNdist lchild_B = bnd_knn[nq_split->child[ANN_LO]].B; // left child bound
      ANNdist rchild_B = bnd_knn[nq_split->child[ANN_HI]].B; // right child bound
      ANNdist b1 = std::min(lchild_B + 2 * (nq_bound.lambda - bounds[nq_split->child[ANN_LO]].lambda),
                            rchild_B + 2 * (nq_bound.lambda - bounds[nq_split->child[ANN_HI]].lambda));
      bound = std::min(b1, N_q_par == NULL ? ANN_DIST_INF : bnd_knn[N_q_par].B);
      // R_INFO(node_labels.at(N_q) << ": min(max(INF), INF, " << b1 << ", " << (N_q_par == NULL ? ANN_DIST_INF : (*bnd_knn)[N_q_par].B) << ")\n")
    } else {
      ANNkd_leaf* nq_leaf = AS_LEAF(N_q); // get leaf node
      ANNdist max_knn = nq_knn_bnd.maxKNN(nq_leaf->n_pts);
      ANNdist b2 = nq_knn_bnd.min_knn + nq_bound.rho + nq_bound.lambda;
      bound = std::min(std::min(max_knn, b2),
                       N_q_par == NULL ? ANN_DIST_INF : bnd_knn[N_q_par].B);
      // R_INFO(node_labels.at(N_q) << ": min(" << max_knn << ", " << b2 << ", INF, " << (N_q_par == NULL ? ANN_DIST_INF : (*bnd_knn)[N_q_par].B) << ")\n")
    }

    // Final bound (lower)
    nq_knn_bnd.B = bound;

    // R_INFO(node_labels.at(N_q) << ": final Bound == " << nq_knn_bnd.B << (IS_LEAF(N_q) ? " (leaf)\n" : "(split)\n"))
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
    // if (N_q != N_r){
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
  }
  return;
}

// Regular (non-pruning) dual tree DFS traversal
// This shouldn't actually be used in practice, but it's a good function to have to benchmark
// how well the dual traversal performs with and without pruning enabled
void DualTreeKNN::DFS(ANNkd_node* N_q, ANNkd_node* N_r){  INC_TRAVERSAL(1)
  // KD Trees only store points in the leaves; Base case only needed if comparing two leaves
  if (IS_LEAF(N_q) && IS_LEAF(N_r)){ ANN_LEAF(2) // Update leaf counter
    BaseCase(N_q, N_r); // Pass nodes as well to keep track of min_knn
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

ANNdist DualTreeKNN::Score(ANNkd_node* N_q, ANNkd_node* N_r) {
  if (N_q == N_r) return 0;
  ANNdist min_dist_qr = min_dist(N_q, N_r); // minimum distance between two bounding rectangles
  ANNdist best_bound = B(N_q); // "best" lower bound
  R_INFO("Min dist. Q (" << node_labels.at(N_q) << ") <--> R (" << node_labels.at(N_r) << "): " << min_dist_qr)
  R_INFO(" (prune if bigger than ==> B(" << node_labels.at(N_q) << "): " << best_bound << " )\n")// : " << max_knn_B(N_q) << ", b2: "
  if (min_dist_qr > best_bound){
    return ANN_DIST_INF; // Prune this branch
  }
  return min_dist_qr; // Recurse into this branch
};

// TODO: rewrite second loop of base case to reset j to track i in the case the nodes are equivalent
ANNdist DualTreeKNN::BaseCaseIdentity(ANNkd_node* N_q, ANNkd_node* N_r){
  ANNdist dist;			// distance to data point
  int q_idx, r_idx; // Indices ot query and reference points
  ANNdist min_dist_q, min_dist_r; // best distances of query and reference points
  ANNkd_leaf* N_q_leaf = AS_LEAF(N_q), *N_r_leaf = AS_LEAF(N_r); // type-checking done previously

  for (int q_i = 0; q_i < N_q_leaf->n_pts; ++q_i){
    for (int r_i = 0; r_i < N_r_leaf->n_pts; ++r_i){
      q_idx = N_q_leaf->bkt[q_i], r_idx = N_r_leaf->bkt[r_i];

      // Compute Base case, saving knn ids and distances along the way
      if (!hasBeenChecked(q_idx, r_idx)) { ANN_PTS(2) // Has this pair been considered before?
        min_dist_q = (knn.at(q_idx))->max_key(); // k-th smallest distance so far (query)
        min_dist_r = (knn.at(r_idx))->max_key(); // k-th smallest distance so far (reference)
        dist = q_idx == r_idx ? 0 : computeDistance(q_idx, r_idx, min_dist_q, min_dist_r);

        R_INFO("min_dist_q: " << min_dist_r << ", min_dist_r: " << min_dist_r << std::endl;)
        R_INFO("new dist: " << dist << " (" << q_idx << ", " << r_idx << ")" << std::endl;)

        // Update the priority queue and cache
        bool add_query_known = false, add_ref_known = false;
        if (dist < min_dist_q){
          add_query_known = (knn.at(q_idx))->insert(dist, r_idx);
          R_INFO("Inserting " << r_idx << " into " << q_idx << " knn list\n")
          R_INFO("Current list (" << q_idx << "): " <<
            (knn.at(q_idx))->ith_smallest_info(0) <<
            (knn.at(q_idx))->ith_smallest_info(1) << "\n")
        }
        if (dist < min_dist_r && q_idx != r_idx){
          add_ref_known = (knn.at(r_idx))->insert(dist, q_idx);
          R_INFO("Inserting " << q_idx << " into " << r_idx << " knn list\n")
          R_INFO("Current list (" << r_idx << "): " <<
            (knn.at(r_idx))->ith_smallest_info(0) <<
            (knn.at(r_idx))->ith_smallest_info(1) << "\n")
        }

        // Refresh the kth smallest again and update bounds
        min_dist_q = (knn.at(q_idx))->max_key(); // k-th smallest distance so far (query)
        min_dist_r = (knn.at(r_idx))->max_key(); // k-th smallest distance so far (reference)
        R_INFO("(new) min_dist_q: " << min_dist_r << ", min_dist_r: " << min_dist_r << std::endl;)
        updateBounds(dist, N_q_leaf, N_r_leaf, min_dist_q, min_dist_r, add_query_known, add_ref_known); // Update KNN bound info.
      } // if(!hasBeenChecked(q_idx, r_idx))
    }
  }
  return dist;
}

ANNdist DualTreeKNN::BaseCaseNonIdentity(ANNkd_node* N_q, ANNkd_node* N_r){
  // TODO
  return 0;
} // end BaseCaseIDC







