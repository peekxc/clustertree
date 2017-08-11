#ifndef DT_KNN_H
#define DT_KNN_H

#include <Rcpp.h>
using namespace Rcpp;

#include <DT/dt.h> // Dual tree class definitions
#include <ANN/structures/pr_queue_k.h> // KNN needs priority queue to record ANNidx and ANNdist records
#include <utilities.h> // R_info, profiling macros, etc.
#include "BoundKNN.h" // Extra bounds specifically for KNN searches

// Useful globals to refer to in the recursion
extern ANNkd_node* N_q_par, *N_r_par;

#ifdef ANN_PERF
extern unsigned int n_traversals;
#define INC_TRAVERSAL(x) n_traversals++;
#define RESET_TRAVERSAL() n_traversals = 0;
#else
#define INC_TRAVERSAL(x)
#define RESET_TRAVERSAL()
#endif

class DualTreeKNN : public DualTree {
protected:
  int k;
  bool knn_identity; // are reference and query sets the same?
  std::unordered_map<ANNidx, ANNmin_k*>* knn;  // Map between point index and kNN priority queue
  std::unordered_map<ANNkd_node*, BoundKNN& >* bnd_knn;
  std::vector<ANNdist>* D; // array of cached knn distances for every point
public:

  // Main constructors
  DualTreeKNN(const bool prune, const int dim); // default constructor
  virtual void setup(ANNkd_tree* kd_treeQ, ANNkd_tree* kd_treeR);
  //~DualTreeKNN();

  // New methods for the derived class
  void KNN(int k, NumericMatrix& dists, IntegerMatrix& ids);
  ANNdist max_knn(ANNkd_node* N_q);

  // Base-class functions replaced for KNN
  // virtual ANNdist min_dist(ANNkd_node* N_i, ANNkd_node* N_j); // use tracked min_knn to recursively compute
  virtual ANNdist B(ANNkd_node* N_q);
  ANNdist max_knn_B(ANNkd_node* N_q);
  virtual void pDFS(ANNkd_node* N_q, ANNkd_node* N_r);
  virtual void DFS(ANNkd_node* N_q, ANNkd_node* N_r);

  // Update various KNN-related bounds
  ANNdist updateBounds(ANNdist new_dist, ANNkd_node* N_q_leaf, ANNkd_node* N_r_leaf, int q_idx, int r_idx);

  // Score function: try to inline if possible
  virtual inline ANNdist Score(ANNkd_node* N_q, ANNkd_node* N_r) {
    if (N_q == N_r) return 0;
    ANNdist min_dist_qr = min_dist(N_q, N_r); // minimum distance between two bounding rectangles
    //ANNdist max_knn_dist = ANN_DIST_INF; // max_knn_B(N_q); // max descendent knn distance of N_q
    ANNdist best_bound = B(N_q); // "best" lower bound
    R_INFO("Min dist. Q (" << node_labels.at(N_q) << ") <--> R (" << node_labels.at(N_r) << "): " << min_dist_qr)
    // R_INFO(" (prune if bigger than ==> b1: " << max_knn_dist << ", b2: " << best_bound << " )\n")

    // ANNdist bound_nq = B(N_q); // Bound type 1 (TODO: doesnt work right now!)
    if (min_dist_qr > best_bound){
      return ANN_DIST_INF; // Prune this branch
    }
    return min_dist_qr; // Recurse into this branch
  };

  // Simple Base case function: try to inline if possible
  virtual inline ANNdist BaseCase(ANNpoint p_q, ANNpoint p_r, const int q_idx, const int r_idx, ANNkd_node* N_q) {
    // Compute distance and retrieve KNN pr queue
    ANNdist dist = annDist(d, qtree->thePoints()[q_idx], rtree->thePoints()[r_idx]);// (ANNdist) ANN_SUM(box_dist, ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));
    ANNmin_k& q_knn = *knn->at(q_idx); // Retrieve query point's knn priority queue

    // Is the computed distance less than the current k-nearest neighbor?
    if (dist < q_knn.max_key()){ q_knn.insert(dist, r_idx); }

    // If the query == reference sets, set the other nearest neighbor as well!
    if (knn_identity && q_idx != r_idx && dist < knn->at(r_idx)->max_key()){
      knn->at(r_idx)->insert(dist, q_idx);
    }
    return dist;
  };

  // Note: *** Minkowski metrics only *** ( to test? )
  inline void BaseCaseNonIdentity(ANNkd_leaf* N_q_leaf, ANNkd_leaf* N_r_leaf){

    ANNdist dist;				// distance to data point
    ANNcoord* pp;				// data coordinate pointer
    ANNcoord* qq;				// query coordinate pointer
    ANNdist min_dist;			// distance to k-th closest point
    ANNcoord t;
    int d_i;

    for (int q_i = 0, q_idx = N_q_leaf->bkt[q_i]; q_i < N_q_leaf->n_pts; ++q_i, q_idx = N_q_leaf->bkt[q_i]){

      ANNmin_k& q_knn = *knn->at(q_idx); // Retrieve query point's knn priority queue
      min_dist = q_knn.max_key(); // k-th smallest distance so far

      for (int r_i = 0; r_i < N_r_leaf->n_pts; ++r_i){
        const int r_idx = N_r_leaf->bkt[r_i];

        // Compute Base case, saving knn ids and distances along the way
        if (!hasBeenChecked(q_idx, r_idx)) { // Has this pair been considered before?
          R_INFO("Calling base case for: q = " << q_idx << ", r = " << r_idx << ")\n")

          // Update number of points known with non-inf knn distances
          (*bnd_knn)[N_q_leaf].knn_known++;

          pp = rtree->thePoints()[r_idx];			// first coord of next data point
          qq = qtree->thePoints()[q_idx];			// first coord of query point
          dist = 0;

          for(d_i = 0; d_i < d; d_i++) {
            ANN_COORD(1)				// one more coordinate hit
            ANN_FLOP(4)					// increment floating ops

            t = *(qq++) - *(pp++);		// compute length and adv coordinate
            // exceeds dist to k-th smallest?
            if( (dist = ANN_SUM(dist, ANN_POW(t))) > min_dist) {
              break;
            }
          }

          if (d_i >= d &&					// among the k best?
              (ANN_ALLOW_SELF_MATCH || dist!=0)) { // and no self-match problem
            // add it to the list
            q_knn.insert(dist, r_idx);
            min_dist = q_knn.max_key();
          }
        }

        // Extra step #1: Update the smallest knn distance for this query node (since we're here)
        ANNmin_k& query_knn = *(*knn)[q_idx], &ref_knn = *(*knn)[r_idx];
        ANNdist knn_dist = N_q_leaf == N_r_leaf ? std::min(ref_knn.max_key(), query_knn.max_key()) : query_knn.max_key();
        if (knn_dist < (*bnd_knn)[N_q_leaf].min_knn){ (*bnd_knn)[N_q_leaf].min_knn = knn_dist; }
        // ANNdist min_knn = std::min(min_dist, N_q_leaf == N_r_leaf ? std::min((*knn)[r_idx]->max_key(), ) : (*bnd_knn)[N_q_leaf].min_knn)
        // (*bnd_knn)[N_q_leaf].min_knn = std::min(min_dist, (*bnd_knn)[N_q_leaf].min_knn);

        // Extra step #2: Check the query pt's and, if applicable, the reference pts knn distance as well to see if
        // the maximum non-inf knn distance should be updated for the current query node. Since we're
        // computing distances in betwen r_idx and q_idx, these keys should be guarenteed to be finite, and there's
        // no need to check how many KNN distances have already been computed.
        BoundKNN& qbnd = bnd_knn->at(N_q_leaf);
        if ((*knn)[q_idx]->max_key() > qbnd.max_real_knn){ qbnd.max_real_knn = (*knn)[q_idx]->max_key(); }
        if (knn_identity && q_idx != r_idx && knn->at(r_idx)->max_key() > qbnd.max_real_knn){
          qbnd.max_real_knn = knn->at(r_idx)->max_key();
        }
        ANN_PTS(N_r_leaf->n_pts) // increment point counter
      }
      ANN_PTS(N_q_leaf->n_pts) // increment point counter
    }
  } // end BaseCaseIDC

  inline ANNdist BaseCaseIdentity(ANNkd_leaf* N_q_leaf, ANNkd_leaf* N_r_leaf){

    ANNdist dist;				// distance to data point
    ANNdist min_dist_q, min_dist_r;			// distance to k-th closest point
    ANNcoord* pp; // first coord of reference point
    ANNcoord* qq; // first coord of query point
    ANNcoord t;
    int d_i;

    // const bool SAME_LEAF = N_q_leaf == N_r_leaf;
    for (int q_i = 0, q_idx = N_q_leaf->bkt[q_i]; q_i < N_q_leaf->n_pts; ++q_i, q_idx = N_q_leaf->bkt[q_i]){

      ANNmin_k& q_knn = (*knn->at(q_idx)); // Retrieve query point's knn priority queue
      min_dist_q = q_knn.max_key(); // k-th smallest distance so far

      for (int r_i = 0, r_idx = N_r_leaf->bkt[r_i]; r_i < N_r_leaf->n_pts; ++r_i, r_idx = N_r_leaf->bkt[r_i]){
        // Compute Base case, saving knn ids and distances along the way
        if (!hasBeenChecked(q_idx, r_idx)) { ANN_PTS(2) // Has this pair been considered before?

          ANNmin_k& r_knn = (*knn->at(r_idx)); // Retrieve ref point's knn priority queue
          min_dist_r = r_knn.max_key(); // k-th smallest distance so far (reference)

          qq = qtree->thePoints()[q_idx];     // first coord of query point
          pp = rtree->thePoints()[r_idx];			// first coord of reference point
          dist = 0;

          // Incrementally compute distance. If at any dimension the distance exceeds the kth smallest
          // distance so far, continue on to the next reference point
          for(d_i = 0; d_i < d; d_i++) {
            ANN_COORD(1)				// one more coordinate hit
            ANN_FLOP(4)					// increment floating ops

            t = *(qq++) - *(pp++);		// compute length and adv coordinate
            // exceeds dist to k-th smallest?
            dist = ANN_SUM(dist, ANN_POW(t));
            if(dist > min_dist_q && dist > min_dist_r) { // check both since workign with tree from identical data set
              break;
            }
          }

          const bool valid_dist = d_i >= d  && (ANN_ALLOW_SELF_MATCH || dist!=0); // ensure is valid distance

          // Update query point KNN distance if better
          if (valid_dist && dist < min_dist_q) {
            (*bnd_knn)[N_q_leaf].knn_known += q_knn.insert(dist, r_idx); // Update number of points known with non-inf knn distances
            min_dist_q = q_knn.max_key(); // k-th smallest distance of query point so far
          }

          // Also update reference point index if better as well, since trees are identical
          if (valid_dist && q_idx != r_idx && dist < min_dist_r){
            (*bnd_knn)[N_r_leaf].knn_known += r_knn.insert(dist, q_idx);
            min_dist_r = r_knn.max_key(); // k-th smallest distance of reference point so far
          }

          // Get KNN Bounds
          BoundKNN& qbnd = bnd_knn->at(N_q_leaf), &rbnd = bnd_knn->at(N_r_leaf);

          // Update minimum kth nearest neighbor distances
          qbnd.min_knn = std::min(qbnd.min_knn, min_dist_q); // Update minimum kth nearest distance on query node
          rbnd.min_knn = std::min(rbnd.min_knn, min_dist_r); // Update minimum kth nearest distance on reference node

          // Update maximum (non-inf) distances as well
          qbnd.max_real_knn = std::max(qbnd.max_real_knn, min_dist_q == ANN_DIST_INF ? 0 : min_dist_q);
          rbnd.max_real_knn = std::max(rbnd.max_real_knn, min_dist_r == ANN_DIST_INF ? 0 : min_dist_r);

          // Informational output
          R_INFO("Q (" << node_labels.at(N_q_leaf) << ") max (non-inf) knn: " << qbnd.max_real_knn
                       << " (" << qbnd.knn_known << "/" << N_q_leaf->n_pts << " known)"
                       << " max_knn: " << qbnd.maxKNN(N_q_leaf->n_pts)
                       << " min_knn: " << qbnd.min_knn << "\n")
          R_INFO("R (" << node_labels.at(N_r_leaf) << ") max (non-inf) knn: " << rbnd.max_real_knn
                       << " (" << rbnd.knn_known << "/" << N_r_leaf->n_pts << " known)"
                       << " max_knn: " << rbnd.maxKNN(N_r_leaf->n_pts)
                       << " min_knn: " << rbnd.min_knn << "\n")
        } // if(!hasBeenChecked(q_idx, r_idx))
      }
    }
    return dist;
  } // end BaseCaseIDC
};

#endif