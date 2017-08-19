#include <Rcpp.h>
using namespace Rcpp;

#include "dtb.h"


DualTreeBoruvka::DualTreeBoruvka(const bool prune, const int dim, const int n, Metric& m)
  : DualTreeKNN(prune, dim, m), CC(UnionFind(n)) {
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

inline ANNdist DualTreeBoruvka::BaseCaseIdentity(ANNkd_node* N_q, ANNkd_node* N_r) {
  ANNdist dist;				// distance to data point
  int q_idx, r_idx; // Indices ot query and reference points
  ANNdist min_dist_q, min_dist_r; // best distances of query and reference points
  ANNkd_leaf* N_q_leaf = AS_LEAF(N_q), *N_r_leaf = AS_LEAF(N_r); // type-checking done previously
  int comp_index; // the component index of the query point

  bool all_same_comp = true; // All points in the same component? Assume positive case
  for (int q_i = 0; q_i < N_q_leaf->n_pts; ++q_i){
    for (int r_i = 0; r_i < N_r_leaf->n_pts; ++r_i){
      int q_idx = N_q_leaf->bkt[q_i], r_idx = N_r_leaf->bkt[r_i];

      // Compute Base case, saving knn ids and distances along the way
      if (!hasBeenChecked(q_idx, r_idx)) { ANN_PTS(2) // Has this pair been considered before?
        min_dist_q = (*knn->at(q_idx)).max_key(); // k-th smallest distance so far (query)
        min_dist_r = (*knn->at(r_idx)).max_key(); // k-th smallest distance so far (reference)
        dist = q_idx == r_idx ? 0 : computeDistance(q_idx, r_idx, min_dist_q, min_dist_r);
        R_PRINTF("q: %d <--> r: %d = %f", q_idx, r_idx, dist)
        R_INFO(", qknn dist: " << min_dist_q << ", rknn dist: " << min_dist_r << "\n")

        // Update the knn priority queue
        bool add_query_known = false, add_ref_known = false;
        if (dist < min_dist_q){
          R_PRINTF("Adding r_idx: %d to q_idx: %d knn queue\n", r_idx, q_idx)
          add_query_known = (*knn->at(q_idx)).insert(dist, r_idx); }
        if (dist < min_dist_r && q_idx != r_idx){
          R_PRINTF("Adding q_idx: %d to r_idx: %d knn queue\n", q_idx, r_idx)
          add_ref_known = (*knn->at(r_idx)).insert(dist, q_idx);
        }

        // Update the kth smallest distances
        min_dist_q = (*knn->at(q_idx)).max_key(); // updated k-th smallest distance so far (query)
        min_dist_r = (*knn->at(r_idx)).max_key(); // updated k-th smallest distance so far (reference)
        updateBounds(dist, N_q_leaf, N_r_leaf, min_dist_q, min_dist_r, add_query_known, add_ref_known); // Upda

        // Main DTB step: update the shortest edge to the nearest component
        if (q_idx != r_idx){
          const int q_comp = CC.Find(q_idx), r_comp = CC.Find(r_idx);
          if (q_comp != r_comp && dist < D[q_comp]){
            N.at(q_comp).from = q_idx;
            N.at(q_comp).to = r_idx;
            R_INFO("Shortest edge: " << q_idx << " <--> " << r_idx << " " << dist << " (old dist: " << D[q_comp] << ")\n")
            D[q_comp] = dist;
          }
          if (q_comp != r_comp && dist < D[r_comp]){
            N.at(r_comp).from = r_idx;
            N.at(r_comp).to = q_idx;
            R_INFO("Shortest edge: " << r_idx << " <--> " << q_idx << " " << dist << " (old dist: " << D[r_comp] << ")\n")
            D[r_comp] = dist;
          }
        }
      } // if(!hasBeenChecked(q_idx, r_idx))
      // Update whether all points in the node are all part of the same component
      // If so, store the CC index. Else store -1.
      const int q_comp = CC.Find(q_idx), r_comp = CC.Find(r_idx);
      all_same_comp = all_same_comp & (q_comp == r_comp); // Re-check w/ Find incase changed
      ALL_CC_SAME[N_q_leaf] = all_same_comp ? q_comp : -1;
      ALL_CC_SAME[N_r_leaf] = all_same_comp ? q_comp : -1;
    }
  }
  return dist;
}

// Score function for DTB. Uses k=1 KNN type bound information.
inline ANNdist DualTreeBoruvka::Score(ANNkd_node* N_q, ANNkd_node* N_r) {
  // if (N_q == N_r) return 0;
  ANNdist min_dist_qr = min_dist(N_q, N_r); // minimum distance between two bounding rectangles
  ANNdist best_bound = B(N_q); // "best" lower bound
  if (min_dist_qr < best_bound){
    R_INFO("CC: " << node_labels.at(N_q) << " == " << ALL_CC_SAME[N_q] << ", ")
    R_INFO("CC: " << node_labels.at(N_r) << " == " << ALL_CC_SAME[N_r] << "\n")
    if (ALL_CC_SAME[N_q] == ALL_CC_SAME[N_r] && ALL_CC_SAME[N_q] >= 0 && ALL_CC_SAME[N_r] >= 0){
      return ANN_DIST_INF; // Prune this branch: never check branches if every point is already in the same comp.
    }
    return (min_dist_qr); // Check for pruning
  }
  return ANN_DIST_INF; // Prune this branch
}

// void DualTreeBoruvka::pDFS(ANNkd_node* N_q, ANNkd_node* N_r){
//   // KD Trees only store points in the leaves; Base case only needed if comparing two leaves
//   if (IS_LEAF(N_q) && IS_LEAF(N_r)) {
//     ANN_LEAF(2) // Update leaf counter
//     BaseCaseIdentity(AS_LEAF(N_q), AS_LEAF(N_r));
//     return; // If at the base case, don't recurse!
//   }
//
//   // --------------- Begin prioritized recursive calls ---------------
//   if (IS_LEAF(N_q) && !IS_LEAF(N_r)){ ANN_LEAF(1); ANN_SPL(1)
//
//     N_r_par = (ANNkd_node*) N_r; // set current reference node as parent, keep same query parent
//     const ANNkd_split* nr_spl = AS_SPLIT(N_r);
//
//     // Sort and visit closer branches
//     NodeScore scores[2] = { NodeScore(N_q, nr_spl->child[ANN_LO], Score(N_q, nr_spl->child[ANN_LO])),
//                             NodeScore(N_q, nr_spl->child[ANN_HI], Score(N_q, nr_spl->child[ANN_HI])) };
//     sort2_sn_stable(scores); // Sort by score using sorting network
//     if (scores[0].score == ANN_DIST_INF) return; else pDFS(scores[0].lhs, scores[0].rhs);
//     if (scores[1].score == ANN_DIST_INF) return; else pDFS(scores[1].lhs, scores[1].rhs);
//
//   } else if (!IS_LEAF(N_q) && IS_LEAF(N_r)){ ANN_SPL(1); ANN_LEAF(1)
//     N_q_par = N_q; // set current query node as parent, keep same reference parent
//     ANNkd_split *const nq_spl = AS_SPLIT(N_q);
//
//     // Sort closer branches
//     NodeScore scores[2] = { NodeScore(nq_spl->child[ANN_LO], N_r, Score(nq_spl->child[ANN_LO], N_r)),
//                             NodeScore(nq_spl->child[ANN_HI], N_r, Score(nq_spl->child[ANN_HI], N_r)) };
//     sort2_sn_stable(scores); // Sort by score using sorting network
//
//     // Visit closer branches first, return at first sign of infinity
//     if (scores[0].score == ANN_DIST_INF) return; else pDFS(scores[0].lhs, scores[0].rhs);
//     if (scores[1].score == ANN_DIST_INF) return; else pDFS(scores[1].lhs, scores[1].rhs);
//
//   } else { ANN_SPL(2)
//     N_q_par = N_q, N_r_par = N_r; // set parents as the current nodes
//     ANNkd_split *const nq_spl = AS_SPLIT(N_q), *const nr_spl = AS_SPLIT(N_r);
//     if (N_q != N_r){
//       NodeScore scores[4] = { NodeScore(nq_spl->child[ANN_LO], nr_spl->child[ANN_LO], Score(nq_spl->child[ANN_LO], nr_spl->child[ANN_LO])),
//                               NodeScore(nq_spl->child[ANN_LO], nr_spl->child[ANN_HI], Score(nq_spl->child[ANN_LO], nr_spl->child[ANN_HI])),
//                               NodeScore(nq_spl->child[ANN_HI], nr_spl->child[ANN_LO], Score(nq_spl->child[ANN_HI], nr_spl->child[ANN_LO])),
//                               NodeScore(nq_spl->child[ANN_HI], nr_spl->child[ANN_HI], Score(nq_spl->child[ANN_HI], nr_spl->child[ANN_HI])) };
//       sort4_sn_stable(scores); // Sort by score using sorting network
//
//       // Visit closer branches first, return at first sign of infinity
//       if (scores[0].score == ANN_DIST_INF) return; else pDFS(scores[0].lhs, scores[0].rhs);
//       if (scores[1].score == ANN_DIST_INF) return; else pDFS(scores[1].lhs, scores[1].rhs);
//       if (scores[2].score == ANN_DIST_INF) return; else pDFS(scores[2].lhs, scores[2].rhs);
//       if (scores[3].score == ANN_DIST_INF) return; else pDFS(scores[3].lhs, scores[3].rhs);
//     } else {
//       // If the nodes are identical, it is unnecessary to score and recurse into both children
//       NodeScore scores[3] = { NodeScore(nq_spl->child[ANN_LO], nr_spl->child[ANN_LO], Score(nq_spl->child[ANN_LO], nr_spl->child[ANN_LO])),
//                               NodeScore(nq_spl->child[ANN_LO], nr_spl->child[ANN_HI], Score(nq_spl->child[ANN_LO], nr_spl->child[ANN_HI])),
//                               NodeScore(nq_spl->child[ANN_HI], nr_spl->child[ANN_HI], Score(nq_spl->child[ANN_HI], nr_spl->child[ANN_HI])) };
//       sort3_sn_stable(scores); // Sort by score using sorting network
//
//       // Visit closer branches first, return at first sign of infinity
//       if (scores[0].score == ANN_DIST_INF) return; else pDFS(scores[0].lhs, scores[0].rhs);
//       if (scores[1].score == ANN_DIST_INF) return; else pDFS(scores[1].lhs, scores[1].rhs);
//       if (scores[2].score == ANN_DIST_INF) return; else pDFS(scores[2].lhs, scores[2].rhs);
//     }
//   }
//   return;
// }
//

List DualTreeBoruvka::DTB(NumericMatrix& x){

  // Necessary data structures for the DTB
  const int n = x.nrow();
  N = std::vector<edgeFT>(n); // one edge per component in the beginning
  D = std::vector<ANNdist>(n, ANN_DIST_INF); // one distance per component

  // Initialize knn structure with k = 2 (1 neighbor outside of identity)
  knn = new std::unordered_map<ANNidx, ANNmin_k*>();
  knn->reserve(n); // reserve enough buckets to hold at least n_pts items
  for (int i = 0; i < n; ++i) {
    knn->insert(std::pair<ANNidx, ANNmin_k*>(qtree->pidx[i], new ANNmin_k(2)));
  }

  // Initialize the components as all false
  ALL_CC_SAME = std::unordered_map<ANNkd_node*, int>();
  ALL_CC_SAME.reserve(n);
  if (use_pruning){
    for (std::unordered_map<ANNkd_node*, const Bound& >::iterator bnd = bounds->begin(); bnd != bounds->end(); ++bnd){
      ALL_CC_SAME.insert(std::make_pair(bnd->first, -1));
    }
  }

  // If pruning is enabled, use that one. Otherwise use regular DFS w/o computing extra bounds.
  NumericMatrix mst_el = NumericMatrix(n - 1, 3);
  if (use_pruning) {
    int c_i = 0;
    while(!fully_connected()){
    //for (int j = 0; j < 1; ++j){
      Rcpp::checkUserInterrupt();

      // Run pruning traversal
      pDFS(rtree->root, qtree->root);

      // Sort the distances
      NumericVector nn_dist = wrap(D);
      IntegerVector index = order_(nn_dist) - 1;

      // Merge components based on nearest edge
      for (IntegerVector::iterator i = index.begin(); i != index.end(); ++i){
        const int from = N.at(*i).from, to = N.at(*i).to;
        if (from != to && D.at(*i) != ANN_DIST_INF && CC.Find(from) != CC.Find(to)){
          mst_el(c_i++, _) = NumericVector::create(from, to, sqrt(D.at(*i)));
          CC.Union(from, to);
        }
      }

      // Reset distances and bounds for next traversal
      std::fill(D.begin(), D.end(), ANN_DIST_INF); // reset shortest edge distance array
      resetBounds(); // Reset KNN bounds
      resetKNN(); // Reset K nearest neighbors priority queue
      resetBaseCases(); // Reset base cases: distances need to be recomputed

      // Informative output
      IntegerVector components = CC.getCC();
      R_INFO("Current CCs: ")
      for (int j = 0; j < n; ++j){ R_INFO(CC.Find(j)); }
      R_INFO("\n")
    }
  } // use_pruning

  IntegerMatrix bc = getBaseCases();
  return List::create(_["mst"] = mst_el, _["base_cases"] = wrap(bc));
}
