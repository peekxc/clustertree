#ifndef DT_KNN_H
#define DT_KNN_H

#include <Rcpp.h>
using namespace Rcpp;

#include <DT/dt.h> // Dual tree class definitions
#include <ANN/structures/pr_queue_k.h> // KNN needs priority queue to record ANNidx and ANNdist records
#include <utilities.h> // R_info, profiling macros, etc.
#include "BoundKNN.h" // Extra bounds specifically for KNN searches

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
  ANNkd_node* N_q_par, *N_r_par;
  int k;
  bool knn_identity; // are reference and query sets the same?
  std::unordered_map<ANNidx, ANNmin_k*>* knn;  // Map between point index and kNN priority queue
  std::unordered_map<ANNkd_node*, BoundKNN& >* bnd_knn;
  std::vector<ANNdist>* D; // array of cached knn distances for every point
public:

  // Main constructors
  DualTreeKNN(const bool prune, const int dim, Metric* m = NULL); // default constructor
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
  ANNdist updateBounds(ANNdist new_dist, ANNkd_node* N_q_leaf, ANNkd_node* N_r_leaf,
                       ANNdist min_dist_q, ANNdist min_dist_r,
                       bool add_known_query, bool add_known_ref);

  // Score function: try to inline if possible
  virtual inline ANNdist Score(ANNkd_node* N_q, ANNkd_node* N_r);
  virtual inline ANNdist BaseCaseIdentity(ANNkd_node* N_q, ANNkd_node* N_r);
  virtual inline ANNdist BaseCaseNonIdentity(ANNkd_node* N_q, ANNkd_node* N_r);

  // inline ANNdist BaseCaseIdentity(ANNkd_leaf* N_q_leaf, ANNkd_leaf* N_r_leaf){
  //
  //   ANNdist dist;				// distance to data point
  //   ANNdist min_dist_q, min_dist_r;			// distance to k-th closest point
  //   ANNcoord* pp; // first coord of reference point
  //   ANNcoord* qq; // first coord of query point
  //   ANNcoord t;
  //   int d_i;
  //
  //   // const bool SAME_LEAF = N_q_leaf == N_r_leaf;
  //   for (int q_i = 0, q_idx = N_q_leaf->bkt[q_i]; q_i < N_q_leaf->n_pts; ++q_i, q_idx = N_q_leaf->bkt[q_i]){
  //
  //     ANNmin_k& q_knn = (*knn->at(q_idx)); // Retrieve query point's knn priority queue
  //     min_dist_q = q_knn.max_key(); // k-th smallest distance so far
  //
  //     for (int r_i = 0, r_idx = N_r_leaf->bkt[r_i]; r_i < N_r_leaf->n_pts; ++r_i, r_idx = N_r_leaf->bkt[r_i]){
  //       // Compute Base case, saving knn ids and distances along the way
  //       if (!hasBeenChecked(q_idx, r_idx)) { ANN_PTS(2) // Has this pair been considered before?
  //
  //         ANNmin_k& r_knn = (*knn->at(r_idx)); // Retrieve ref point's knn priority queue
  //         min_dist_r = r_knn.max_key(); // k-th smallest distance so far (reference)
  //
  //         qq = qtree->thePoints()[q_idx];     // first coord of query point
  //         pp = rtree->thePoints()[r_idx];			// first coord of reference point
  //         dist = 0;
  //
  //         // Incrementally compute distance. If at any dimension the distance exceeds the kth smallest
  //         // distance so far, continue on to the next reference point
  //         for(d_i = 0; d_i < d; d_i++) {
  //           ANN_COORD(1)				// one more coordinate hit
  //           ANN_FLOP(4)					// increment floating ops
  //
  //           t = *(qq++) - *(pp++);		// compute length and adv coordinate
  //           // exceeds dist to k-th smallest?
  //           dist = ANN_SUM(dist, ANN_POW(t));
  //           if(dist > min_dist_q && dist > min_dist_r) { // check both since workign with tree from identical data set
  //             break;
  //           }
  //         }
  //
  //         const bool valid_dist = d_i >= d  && (ANN_ALLOW_SELF_MATCH || dist!=0); // ensure is valid distance
  //
  //         // Update query point KNN distance if better
  //         if (valid_dist && dist < min_dist_q) {
  //           (*bnd_knn)[N_q_leaf].knn_known += q_knn.insert(dist, r_idx); // Update number of points known with non-inf knn distances
  //           min_dist_q = q_knn.max_key(); // k-th smallest distance of query point so far
  //         }
  //
  //         // Also update reference point index if better as well, since trees are identical
  //         if (valid_dist && q_idx != r_idx && dist < min_dist_r){
  //           (*bnd_knn)[N_r_leaf].knn_known += r_knn.insert(dist, q_idx);
  //           min_dist_r = r_knn.max_key(); // k-th smallest distance of reference point so far
  //         }
  //
  //         // Get KNN Bounds
  //         BoundKNN& qbnd = bnd_knn->at(N_q_leaf), &rbnd = bnd_knn->at(N_r_leaf);
  //
  //         // Update minimum kth nearest neighbor distances
  //         qbnd.min_knn = std::min(qbnd.min_knn, min_dist_q); // Update minimum kth nearest distance on query node
  //         rbnd.min_knn = std::min(rbnd.min_knn, min_dist_r); // Update minimum kth nearest distance on reference node
  //
  //         // Update maximum (non-inf) distances as well
  //         qbnd.max_real_knn = std::max(qbnd.max_real_knn, min_dist_q == ANN_DIST_INF ? 0 : min_dist_q);
  //         rbnd.max_real_knn = std::max(rbnd.max_real_knn, min_dist_r == ANN_DIST_INF ? 0 : min_dist_r);
  //
  //         // Informational output
  //         R_INFO("Q (" << node_labels.at(N_q_leaf) << ") max (non-inf) knn: " << qbnd.max_real_knn
  //                      << " (" << qbnd.knn_known << "/" << N_q_leaf->n_pts << " known)"
  //                      << " max_knn: " << qbnd.maxKNN(N_q_leaf->n_pts)
  //                      << " min_knn: " << qbnd.min_knn << "\n")
  //         R_INFO("R (" << node_labels.at(N_r_leaf) << ") max (non-inf) knn: " << rbnd.max_real_knn
  //                      << " (" << rbnd.knn_known << "/" << N_r_leaf->n_pts << " known)"
  //                      << " max_knn: " << rbnd.maxKNN(N_r_leaf->n_pts)
  //                      << " min_knn: " << rbnd.min_knn << "\n")
  //       } // if(!hasBeenChecked(q_idx, r_idx))
  //     }
  //   }
  //   return dist;
  // } // end BaseCaseIDC

  };

#endif