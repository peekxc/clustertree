#include <Rcpp.h>
using namespace Rcpp;

#include "dtb.h"

DualTreeBoruvka::DualTreeBoruvka(const bool prune, const int dim, const int n)
  : DualTree(prune, dim), CC(UnionFind(n)) {
  N_q_par = N_r_par = NULL;
  if (prune){
    R_INFO("Allocating KNN bound map\n")
    bnd_knn = new std::unordered_map<ANNkd_node*, BoundKNN& >();
  }
}

void DualTreeBoruvka::setup(ANNkd_tree* kd_treeQ, ANNkd_tree* kd_treeR) {
  setTrees(kd_treeQ, kd_treeR);  // Call parent class setup to assign trees
  if (use_pruning){
    for (std::unordered_map<ANNkd_node*, const Bound& >::iterator bnd = bounds->begin(); bnd != bounds->end(); ++bnd){
      bnd_knn->insert(std::pair<ANNkd_node*, BoundKNN& >(bnd->first, (BoundKNN&) *new BoundKNN()));
    }
  }
}

inline ANNdist DualTreeBoruvka::BaseCaseIdentity(ANNkd_leaf* N_q_leaf, ANNkd_leaf* N_r_leaf) {

  ANNdist dist;				// distance to data point
  int q_idx, r_idx; // Indices ot query and reference points
  ANNdist min_dist_q, min_dist_r; // best distances of query and reference points
  bool all_same_comp = true; // Assume positive case

  for (int q_i = 0; q_i < N_q_leaf->n_pts; ++q_i]){
    for (int r_i = 0; r_i < N_r_leaf->n_pts; ++r_i){
      int q_idx = N_q_leaf->bkt[q_i], r_idx = N_r_leaf->bkt[r_i];
      if (q_idx == r_idx) continue;

      // Compute Base case, saving knn ids and distances along the way
      if (!hasBeenChecked(q_idx, r_idx)) { ANN_PTS(2) // Has this pair been considered before?
        min_dist_q = EL[q_idx].weight; // k-th smallest distance so far (query)
        min_dist_r = EL[r_idx].weight; // k-th smallest distance so far (reference)
        dist = computeDistance(q_idx, r_idx, min_dist_q, min_dist_r);
      }

      // Update KNN bounds
      updateBounds(dist, N_q_leaf, N_r_leaf, q_idx, r_idx);

      // Update edge weight if closer
      if (dist < min_dist_q) { EL[q_idx].weight = dist; }
      if (dist < min_dist_r) { EL[r_idx].weight = dist; }

      // Update Component indexes
      const bool same_comp = CC.Find(q_idx) == CC.Find(r_idx);
      ALL_CC_SAME[N_q_leaf] = all_same_comp & same_comp;
      ALL_CC_SAME[N_r_leaf] = all_same_comp & same_comp;
      } // if(!hasBeenChecked(q_idx, r_idx))
    }
  }
  return dist;
}

// TODO
inline ANNdist DualTreeBoruvka::Score(ANNkd_node* N_q, ANNkd_node* N_r) {
  if (N_q == N_r) return 0;
  ANNdist min_dist_qr = min_dist(N_q, N_r); // minimum distance between two bounding rectangles
  ANNdist best_bound = B(N_q); // "best" lower bound
  if (min_dist_qr < best_bound){
    if (ALL_CC_SAME[N_q] && ALL_CC_SAME[N_r])
    ANNkd_split* nq_spl = AS_SPLIT(N_q);
    ANNkd_split* nr_spl = AS_SPLIT(N_q);
    return ANN_DIST_INF; // Prune this branch
  }
  return min_dist_qr; // Recurse into this branch
}

ANNdist DualTreeBoruvka::B(ANNkd_node* N_q){
  // If B has been computed before, return it, otherwise look at what needs computed
  BoundKNN& nq_knn_bnd = (BoundKNN&) (*bnd_knn)[N_q];
  if (nq_knn_bnd.B != ANN_DIST_INF){ return(nq_knn_bnd.B); }
  else { // Check to see if we can get better bound

    // Retrieve regular bound object
    const Bound& nq_bound = (const Bound&) (*bounds)[N_q];

    // There are 4 bounds to compute, although some apply only to leaf or split nodes
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
    return nq_knn_bnd.B;
  }
}

void DualTreeBoruvka::pDFS(ANNkd_node* N_q, ANNkd_node* N_r){
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

    // Sort and visit closer branches
    NodeScore scores[2] = { NodeScore(N_q, nr_spl->child[ANN_LO], Score(N_q, nr_spl->child[ANN_LO])),
                            NodeScore(N_q, nr_spl->child[ANN_HI], Score(N_q, nr_spl->child[ANN_HI])) };
    sort2_sn_stable(scores); // Sort by score using sorting network
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
      if (scores[0].score == ANN_DIST_INF) return; else pDFS(scores[0].lhs, scores[0].rhs);
      if (scores[1].score == ANN_DIST_INF) return; else pDFS(scores[1].lhs, scores[1].rhs);
      if (scores[2].score == ANN_DIST_INF) return; else pDFS(scores[2].lhs, scores[2].rhs);
    }
  }
  return;
}


NumericMatrix DualTreeBoruvka::DTB(NumericMatrix& x){

  // New edge list
  const int n = x.nrow();
  EL = new double_edge[n]; // edge list
  for (int i = 0; i < n; ++i){ EL[i].weight = ANN_DIST_INF; }
  ALL_CC_SAME = new std::unordered_map<ANNkd_node*, bool>();
  ALL_CC_SAME.reserve(n);

  // Initialize the components as all false
  if (use_pruning){
    const bool init_value = qtree->bkt_size == 1;
    for (std::unordered_map<ANNkd_node*, const Bound& >::iterator bnd = bounds->begin(); bnd != bounds->end(); ++bnd){
      ALL_CC_SAME.insert(std::make_pair(bnd->first, init_value));
    }
  }

  // If pruning is enabled, use that one. Otherwise use regular DFS w/o computing extra bounds.
  pDFS(rtree->root, qtree->root);

  // Copy over to MST
  NumericMatrix mst_el = NumericMatrix(n, 3);
  for (int i = 0; i < n; ++i){
    mst_el(i, _) = NumericVector::create(EL[i].from, EL[i].to, EL[i].weight);
  }
  return mst_el;
}
