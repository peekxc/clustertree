#include <Rcpp.h>
using namespace Rcpp;

#include "dtb.h"


DualTreeBoruvka::DualTreeBoruvka(const NumericMatrix& q_x, Metric& m, NumericMatrix& r_x, List& config)
  : DualTreeKNN(q_x, m, r_x, config), CC(UnionFind(q_x.nrow())) {
  // No initialization to do beyond this
}

// Base case for DTB (when the query and reference sets are identical)
inline ANNdist DualTreeBoruvka::BaseCaseIdentity(ANNkd_node* N_q, ANNkd_node* N_r) {
  ANNdist dist;				// distance to data point
  int q_idx, r_idx; // Indices ot query and reference points
  ANNdist min_dist_q, min_dist_r; // best distances of query and reference points
  ANNkd_leaf* N_q_leaf = AS_LEAF(N_q), *N_r_leaf = AS_LEAF(N_r); // type-checking done previously
  //int comp_index; // the component index of the query point

  bool all_same_comp = true; // All points in the same component? Assume positive case
  for (int q_i = 0; q_i < N_q_leaf->n_pts; ++q_i){
    for (int r_i = 0; r_i < N_r_leaf->n_pts; ++r_i){
      q_idx = N_q_leaf->bkt[q_i], r_idx = N_r_leaf->bkt[r_i];

      // Compute Base case, saving knn ids and distances along the way
      if (!hasBeenChecked(q_idx, r_idx)) { ANN_PTS(2) // Has this pair been considered before?
        min_dist_q = (knn.at(q_idx))->max_key(); // k-th smallest distance so far (query)
        min_dist_r = (knn.at(r_idx))->max_key(); // k-th smallest distance so far (reference)
        dist = q_idx == r_idx ? 0 : computeDistance(q_idx, r_idx, min_dist_q, min_dist_r);
        R_PRINTF("q: %d <--> r: %d = %f", q_idx, r_idx, dist)
        R_INFO(", qknn dist: " << min_dist_q << ", rknn dist: " << min_dist_r << "\n")

        // Update the knn priority queue
        bool add_query_known = false, add_ref_known = false;
        if (dist < min_dist_q){
          R_PRINTF("Adding r_idx: %d to q_idx: %d knn queue\n", r_idx, q_idx)
          add_query_known = (knn.at(q_idx))->insert(dist, r_idx); }
        if (dist < min_dist_r && q_idx != r_idx){
          R_PRINTF("Adding q_idx: %d to r_idx: %d knn queue\n", q_idx, r_idx)
          add_ref_known = (knn.at(r_idx))->insert(dist, q_idx);
        }

        // Update the kth smallest distances
        min_dist_q = (knn.at(q_idx))->max_key(); // updated k-th smallest distance so far (query)
        min_dist_r = (knn.at(r_idx))->max_key(); // updated k-th smallest distance so far (reference)
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
ANNdist DualTreeBoruvka::Score(ANNkd_node* N_q, ANNkd_node* N_r) {
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

// Dual Tree Boruvka Main function
List DualTreeBoruvka::DTB(const NumericMatrix& x){

  // Necessary data structures for the DTB
  const int n = x.nrow();
  N = std::vector<edgeFT>(n); // one edge per component in the beginning
  D = std::vector<ANNdist>(n, ANN_DIST_INF); // one distance per component

  // Initialize knn structure with k = 2 (1 neighbor outside of identity)
  knn = std::unordered_map<ANNidx, ANNmin_k*>();
  knn.reserve(n); // reserve enough buckets to hold at least n_pts items
  for (int i = 0; i < n; ++i) {
    knn.insert(std::pair<ANNidx, ANNmin_k*>(qtree->pidx[i], new ANNmin_k(2)));
  }

  // Initialize the components as all false
  ALL_CC_SAME = std::unordered_map<ANNkd_node*, int>();
  ALL_CC_SAME.reserve(n);
  if (use_pruning){
    for (std::unordered_map<ANNkd_node*, const Bound& >::iterator bnd = bounds.begin(); bnd != bounds.end(); ++bnd){
      ALL_CC_SAME.insert(std::make_pair(bnd->first, -1));
    }
  }

  // If pruning is enabled, use that one. Otherwise use regular DFS w/o computing extra bounds.
  IntegerMatrix mst_el = IntegerMatrix(n - 1, 2);
  NumericVector height = NumericVector(n - 1);
  int c_i = 0;
  if (use_pruning) {
    // IntegerVector old_CC, new_CC;
    do {
      Rcpp::checkUserInterrupt();

      // Run pruning traversal
      pDFS(rtree->root, qtree->root);

      // Sort the distances
      NumericVector nn_dist = wrap(D); // D is updated with the k=2 nearest neighbor distances from the boruvka step
      IntegerVector index = order_(nn_dist) - 1; // Prepare to connect the small edges first

      // Merge components based on nearest edge
      for (IntegerVector::iterator i = index.begin(); i != index.end(); ++i){
        const int from = N.at(*i).from, to = N.at(*i).to;
        if (from != to && D.at(*i) != ANN_DIST_INF && CC.Find(from) != CC.Find(to)){
          mst_el(c_i, _) = IntegerVector::create(from, to); // Record the point ids
          height.at(c_i) = D.at(*i); // Record the (unfinalized!) distance
          CC.Union(from, to); // Union the two points
          c_i++;
        }
      }

      // Reset distances and bounds for next traversal
      std::fill(D.begin(), D.end(), ANN_DIST_INF); // reset shortest edge distance array
      resetBounds(); // Reset KNN bounds
      resetKNN(); // Reset K nearest neighbors priority queue
      resetBaseCases(); // Reset base cases: distances need to be recomputed

    // Continue until either fully connected or until the
    } while(c_i != n && !fully_connected());

    // Debugging output
    #ifdef NDEBUG
        R_INFO("Final CC: ");
        CC.printCC();
    #endif

  } // use_pruning
  // IntegerVector valid_idx = Rcpp::Range(0, c_i - 1);(Rcpp::Range(0, c_i - 1), Rcpp::Range(0, 1))


  // Replace invalid-distances with NAs
  for (int i = 0; i < n - 1; ++i){
    if (height[i] == std::numeric_limits<ANNdist>::max()){
      height[i] = NA_REAL;
    }
  }
  m_dist.finalize(height); // apply final metric transformation
  return List::create(_["mst"] = mst_el, _["height"] = height, _["CCs"] = CC.getCC());
}




