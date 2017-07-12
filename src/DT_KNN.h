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
  BoundKNN() : B(-1.0), knn_known(0), min_knn(ANN_DIST_INF), max_real_knn(0) {};
  // ANNdist maxKNN(){}
};

class DualTreeKNN : public DualTree {
protected:
  bool knn_identity; // are reference and query sets the same?
  std::unordered_map<ANNidx, ANNmin_k*>* knn;  // Map between point index and kNN priority queue
  std::unordered_map<ANNkd_node*, BoundKNN& >* bnd_knn;
  // std::vector<ANNidx>* qpts;
public:

  // Main constructors
  DualTreeKNN(const bool prune, const int dim); // default constructor
  virtual void setup(ANNkd_tree* kd_treeQ, ANNkd_tree* kd_treeR);
  //~DualTreeKNN();

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
  virtual inline ANNdist Score(ANNkd_node* N_q, ANNkd_node* N_r) final {
      R_INFO("Scoring: Q = " << N_q << ", " <<  "R = " << N_r << "\n")
      ANNdist min_dist_qr = min_dist(N_q, N_r);
      R_INFO("Min dist. Q <--> R: " << sqrt(min_dist_qr) << "\n")

      // ANNdist bound_nq = B(N_q); // Bound type 1 (TODO: doesnt work right now!)
      ANNdist bound_nq = max_knn_B(N_q); // Bound type 2

      if (min_dist_qr < bound_nq){ // check if the minimum distance between two nodes is small enough to recurse
        return (min_dist_qr);
      }
      R_INFO("Pruning! Q <--> R: " << min_dist_qr << "(> " << bound_nq << ")" <<  "\n")
      return ANN_DIST_INF; // Prune this branch
      // return 0;// This branch must be recursed into
  };

  // Base case function: try to inline if possible
  virtual inline void BaseCase(ANNpoint p_q, ANNpoint p_r, const int q_idx, const int r_idx, ANNkd_node* N_q) final{
    R_PRINTF("dist(%d, %d): ", q_idx, r_idx);
    if (!((bool) ((*BC_check)[std::minmax(q_idx, r_idx)]))){ // Has this pair been considered before?

      // Update number of points known with non-inf knn distances
      if (use_pruning) { (*bnd_knn)[N_q].knn_known += (knn_identity && q_idx != r_idx) ? 2 : 1; }


      ANNdist dist = annDist(d, qtree->pts[q_idx], rtree->pts[r_idx]);// (ANNdist) ANN_SUM(box_dist, ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));
      R_PRINTF("%f \n", dist);

      ANNmin_k& q_knn = *knn->at(q_idx); // Retrieve query point's knn

      // Is the computed distance less than the current k-nearest neighbor?
      if (dist < q_knn.max_key()){
        q_knn.insert(dist, r_idx);
        R_INFO("Storing better neighbor for: " << q_idx << "( id=" << r_idx << ", dist=" << dist << ")\n")
      }

      // If the query == reference sets, set the other nearest neighbor as well!
      if (knn_identity && q_idx != r_idx && dist < knn->at(r_idx)->max_key()){
        knn->at(r_idx)->insert(dist, q_idx);
        R_INFO("Storing better neighbor for: " << r_idx << "( id=" << q_idx << ", dist=" << dist << ")\n")
      }

      (*BC_check)[std::minmax(q_idx, r_idx)] = true; // Equivalent pairs won't be visited
      if (!use_pruning) { return; }
      else {
      // Check the query pt's and, if applicable, the reference pts knn distance as well to see if
      // the maximum non-inf knn distance should be updated for the current query node. Since we're
      // computing distances in betwen r_idx and q_idx, these keys should be guarenteed to not be
      // infinite.
        BoundKNN& qnode_knn_bnd = bnd_knn->at(N_q);
        if (knn_identity && q_idx != r_idx && knn->at(r_idx)->max_key() > qnode_knn_bnd.max_real_knn){
          assert(knn->at(r_idx)->max_key() != ANN_DIST_INF);
          qnode_knn_bnd.max_real_knn = knn->at(r_idx)->max_key();
        }
        if (q_knn.max_key() > qnode_knn_bnd.max_real_knn){
          assert(q_knn.max_key() != ANN_DIST_INF);
          qnode_knn_bnd.max_real_knn = q_knn.max_key();
        }
      }
    }
  };
};