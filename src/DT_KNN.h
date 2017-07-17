#include <Rcpp.h>
using namespace Rcpp;

#include "dt.h"
#include "pr_queue_k.h" // KNN needs priority queue to record ANNidx and ANNdist records
#include <array>

// Pruning DT-DFS only: using regular Bound object for tree bounds, but
// provide additional bounds specific to KNN
struct BoundKNN {
  ANNdist B; // Actual computed bound
  unsigned int knn_known; // the number of descendent (or held) points that have at least 1 non-inf KNN distance
  ANNdist min_knn; // minimum distance to the kth nearest neighbor of the points within the node
  ANNdist max_real_knn; // maximum non-inf knn distance
  BoundKNN() : B(ANN_DIST_INF), knn_known(0), min_knn(ANN_DIST_INF), max_real_knn(0) {};
  // ANNdist maxKNN(){}
};

class DualTreeKNN : public DualTree {
protected:
  bool knn_identity; // are reference and query sets the same?
  std::unordered_map<ANNidx, ANNmin_k*>* knn;  // Map between point index and kNN priority queue
  std::unordered_map<ANNkd_node*, BoundKNN& >* bnd_knn;
  // ScoreMap pr_scores;
  // std::vector<ANNidx>* qpts;
public:

  // Main constructors
  DualTreeKNN(const bool prune, const int dim); // default constructor
  virtual void setup(ANNkd_tree* kd_treeQ, ANNkd_tree* kd_treeR);
  //~DualTreeKNN();

  inline const bool hasBeenChecked(const ANNidx q_idx, const ANNidx r_idx){
    // inserts new entry (w/ value false) if the value didn't exist
    return (const bool) (*BC_check)[std::minmax(q_idx, r_idx)];
  }

  inline void hasBeenChecked(const ANNidx q_idx, const ANNidx r_idx, const bool status){
   (*BC_check)[std::minmax(q_idx, r_idx)] = status;
  }

  // New methods for the derived class
  void KNN(int k, NumericMatrix& dists, IntegerMatrix& ids);
  ANNdist max_knn(ANNkd_node* N_q);

  // Base-class functions replaced for KNN
  virtual ANNdist min_dist(ANNkd_node* N_i, ANNkd_node* N_j); // use tracked min_knn to recursively compute
  virtual ANNdist B(ANNkd_node* N_q);
  ANNdist max_knn_B(ANNkd_node* N_q);
  virtual void pDFS(ANNkd_node* N_q, ANNkd_node* N_r);
  virtual void DFS(ANNkd_node* N_q, ANNkd_node* N_r);

  // Score function: try to inline if possible
  virtual inline ANNdist Score(ANNkd_node* N_q, ANNkd_node* N_r) {
    if (N_q == N_r) return 0;
    // R_INFO("Scoring: Q = " << N_q << ", " <<  "R = " << N_r << "\n")

    // The minimum distance between two nodes is defined as the distance between the centroids of the
    // two smallest spheres fitting inside the boudning box of a given node. Since the bounding boxes are
    // n-dimensional rectangles, this is simply the shortest difference between the bounds of each dimension
    // divided by 2.
    ANNdist min_dist_qr = min_dist(N_q, N_r), max_knn_dist = max_knn_B(N_q);
    R_INFO("Min dist. Q (" << N_q << ") <--> R (" << N_r << "): " << min_dist_qr << " (smaller than? " << max_knn_dist << " )\n")

    // ANNdist bound_nq = B(N_q); // Bound type 1 (TODO: doesnt work right now!)
    // ANNdist bound_nq =(*bnd_knn)[N_q].max_real_knn ; // max_knn_B(N_q); // Bound type 2

    // Check if the minimum distance between two nodes is smaller (or small enough) compared to
    // the maximum descendent KNN distance. If so, this pair needs to be recursed into.
    if (min_dist_qr < max_knn_dist){
      return (min_dist_qr); // This branch must be recursed into
    }

    // If the maximum knn distance of the query node is above the minimum distance between the ref. and
    // query nodes (and, by extension, all of their desc. points), no desc. point of the reference node
    // could possibly be the nearest neighbor of the descendents of the query node, so prune this branch.
    // R_INFO("Pruning! Q <--> R: " << min_dist_qr << "( > " << (*bnd_knn)[N_q].max_real_knn << ")" <<  "\n")
    return ANN_DIST_INF; // Prune this branch
    // return 0;// This branch must be recursed into
  };

  // Base case function: try to inline if possible
  virtual inline void BaseCase(ANNpoint p_q, ANNpoint p_r, const int q_idx, const int r_idx, ANNkd_node* N_q) {

  ANNdist dist;				// distance to data point
  ANNcoord* pp;				// data coordinate pointer
  ANNcoord* qq;				// query coordinate pointer
  ANNdist min_dist;			// distance to k-th closest point
  ANNcoord t;
  int d;

  ANNmin_k& q_knn = *knn->at(q_idx); // Retrieve query point's knn priority queue
  min_dist = q_knn.max_key(); // k-th smallest distance so far

  for (int i = 0; i < n_pts; i++) {	// check points in bucket

    pp = ANNkdPts[bkt[i]];			// first coord of next data point
    qq = ANNkdQ;					// first coord of query point
    dist = 0;

    for(d = 0; d < ANNkdDim; d++) {
      ANN_COORD(1)				// one more coordinate hit
      ANN_FLOP(4)					// increment floating ops

      t = *(qq++) - *(pp++);		// compute length and adv coordinate
      // exceeds dist to k-th smallest?
      if( (dist = ANN_SUM(dist, ANN_POW(t))) > min_dist) {
        break;
      }
    }

    if (d >= ANNkdDim &&					// among the k best?
        (ANN_ALLOW_SELF_MATCH || dist!=0)) { // and no self-match problem
      // add it to the list
      ANNkdPointMK->insert(dist, bkt[i]);
      min_dist = ANNkdPointMK->max_key();
    }
  }

    // // Compute distance and retrieve KNN pr queue
    // ANNdist dist = annDist(d, qtree->pts[q_idx], rtree->pts[r_idx]);// (ANNdist) ANN_SUM(box_dist, ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));
    // ANNmin_k& q_knn = *knn->at(q_idx); // Retrieve query point's knn priority queue
    //
    // // Is the computed distance less than the current k-nearest neighbor?
    // if (dist < q_knn.max_key()){
    //   q_knn.insert(dist, r_idx);
    //   // R_INFO("Storing better neighbor for: " << q_idx << " (id=" << r_idx << ", dist=" << dist << ")\n")
    // }
    //
    // // If the query == reference sets, set the other nearest neighbor as well!
    // if (knn_identity && q_idx != r_idx && dist < knn->at(r_idx)->max_key()){
    //   knn->at(r_idx)->insert(dist, q_idx);
    //   // R_INFO("Storing better neighbor for: " << r_idx << " (id=" << q_idx << ", dist=" << dist << ")\n")
    // }
  };


  // Base case using incremental distance calculations
  // Note: *** Minkowski metrics only ***
  inline void BaseCaseIDC(ANNkd_leaf* N_q_leaf, ANNkd_leaf* N_r_leaf){

    ANNdist dist;				// distance to data point
    ANNcoord* pp;				// data coordinate pointer
    ANNcoord* qq;				// query coordinate pointer
    ANNdist min_dist;			// distance to k-th closest point
    ANNcoord t;
    int d;

    for (int q_i = 0, q_idx = N_q_leaf->bkt[q_i]; q_i < N_q_leaf->n_pts; ++q_i, q_idx = N_q_leaf->bkt[q_i]){

      ANNmin_k& q_knn = *knn->at(q_idx); // Retrieve query point's knn priority queue
      min_dist = q_knn.max_key(); // k-th smallest distance so far

      for (int r_i = 0; r_i < N_r_leaf->n_pts; ++r_i){
        const int r_idx = N_r_leaf->bkt[r_i];

        // Compute Base case, saving knn ids and distances along the way
        if (!hasBeenChecked(q_idx, r_idx)) { // Has this pair been considered before?

          for (int i = 0; i < n_pts; i++) {	// check points in bucket

            pp = ANNkdPts[bkt[i]];			// first coord of next data point
            qq = ANNkdQ;					// first coord of query point
            dist = 0;

            for(d = 0; d < ANNkdDim; d++) {
              ANN_COORD(1)				// one more coordinate hit
              ANN_FLOP(4)					// increment floating ops

              t = *(qq++) - *(pp++);		// compute length and adv coordinate
              // exceeds dist to k-th smallest?
              if( (dist = ANN_SUM(dist, ANN_POW(t))) > min_dist) {
                break;
              }
            }

            if (d >= ANNkdDim &&					// among the k best?
                (ANN_ALLOW_SELF_MATCH || dist!=0)) { // and no self-match problem
              // add it to the list
              ANNkdPointMK->insert(dist, bkt[i]);
              min_dist = ANNkdPointMK->max_key();
            }
          }



          R_INFO("Calling base case for: q = " << q_idx << ", r = " << r_idx << ")\n")

          // Update number of points known with non-inf knn distances
          (*bnd_knn)[N_q_leaf].knn_known += (knn_identity && q_idx != r_idx) ? 2 : 1;

          // Call the base case, updating the KNN distances where necessary.
          // BaseCase(qtree->pts[q_idx], rtree->pts[r_idx], q_idx, r_idx, N_q_leaf); // Pass nodes as well to keep track of min_knn
          hasBeenChecked(q_idx, r_idx, true); // Equivalent pairs won't be visited again

          // Extra step #1: Update the smallest knn distance for this query node (since we're here)
          ANNmin_k& query_knn = *(*knn)[q_idx], &ref_knn = *(*knn)[r_idx];
          ANNdist knn_dist = N_q == N_r ? std::min(ref_knn.max_key(), query_knn.max_key()) : query_knn.max_key();
          if (knn_dist < (*bnd_knn)[N_q].min_knn){ (*bnd_knn)[N_q].min_knn = knn_dist; }

          // Extra step #2: Check the query pt's and, if applicable, the reference pts knn distance as well to see if
          // the maximum non-inf knn distance should be updated for the current query node. Since we're
          // computing distances in betwen r_idx and q_idx, these keys should be guarenteed to be finite, and there's
          // no need to check how many KNN distances have already been computed.
          BoundKNN& qbnd = bnd_knn->at(N_q);
          if ((*knn)[q_idx]->max_key() > qbnd.max_real_knn){ qbnd.max_real_knn = (*knn)[q_idx]->max_key(); }
          if (knn_identity && q_idx != r_idx && knn->at(r_idx)->max_key() > qbnd.max_real_knn){
            qbnd.max_real_knn = knn->at(r_idx)->max_key();
          }
        }
        ANN_PTS(N_r_leaf->n_pts) // increment point counter
      }
      ANN_PTS(N_q_leaf->n_pts) // increment point counter
    }
  }


};

