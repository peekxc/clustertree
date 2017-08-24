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
  std::unordered_map<ANNidx, ANNmin_k*> knn;  // Map between point index and kNN priority queue
  std::unordered_map<ANNkd_node*, BoundKNN& > bnd_knn;
  //std::vector<ANNdist>* D; // array of cached knn distances for every point
public:

  // Main constructors
  DualTreeKNN(const NumericMatrix& q_x, Metric& m, NumericMatrix& r_x = emptyMatrix, List& config = default_params);
  ~DualTreeKNN();

  // New methods for the derived class
  List KNN(int k, bool finalize = true);
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
  void resetBounds();
  void resetKNN();

  // Score function: try to inline if possible
  virtual ANNdist Score(ANNkd_node* N_q, ANNkd_node* N_r) override;
  virtual ANNdist BaseCaseIdentity(ANNkd_node* N_q, ANNkd_node* N_r) override;
  virtual ANNdist BaseCaseNonIdentity(ANNkd_node* N_q, ANNkd_node* N_r) override;

  };

#endif