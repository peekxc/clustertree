#ifndef DUALTREE_KNN_H
#define DUALTREE_KNN_H
// DualTreeKNN.h
#include "DualTreeSearch.h"

template <class METRIC_T>
class DualTreeKNN : public DualTreeSearch<METRIC_T> {
private:
  typedef DualTreeSearch<METRIC_T> dts; // dual tree search parent/base type definition
public:
  const int k;  // how many k-points to find the nearest neighbors of
  std::vector<ANNmin_k*> knn_pq; // knn priority queue for each point
  // struct knn_cache {
  //   std::vector<ANNdist>* dist_ptr; // pointer to cached distances
  //   std::vector<ANNmin_k::mk_node*> knn_ptr; // pointer to mk_node pointers that point to the knn distances
  // };
  //std::unordered_map< std::pair<ANNkd_node*, ANNkd_node*>, knn_cache > cache;
  // std::unordered_map<ANNkd_node*, std::vector<ANNmin_k::mk_node*> > knn_ptr; // pointer to mk_node pointers that point to the knn distances

  DualTreeKNN() = default; // default constructor allows specific template instantiations
  DualTreeKNN(ANNkd_tree* qtree, ANNkd_tree* rtree, const NODE_INFO& qinfo, const NODE_INFO& rinfo, METRIC_T& m, const int k);
  //   : DualTreeSearch<METRIC_T>(m), k(k) {
  //     knn_pq = std::vector<ANNmin_k*>(qtree->n_pts); // priority queue for the closest k points
  //     for (int i = 0; i < qtree->n_pts; ++i){ knn_pq[i] = new ANNmin_k(k); } // initialize pr queue
  //
  //     // Benchmarking / Profiling
  //     benchmarks = std::vector<double>(9, 0.0); // 0,
  //     dts::n_comparisons = 0;
  //     dts::n_pruned = 0;
  // };

  // Update function for split nodes
  void UpdateBounds(ANNkd_split* node, NODE_INFO& ninfo) override;
  void UpdateBounds(ANNkd_leaf* node, NODE_INFO& ninfo) override;

  // KNN Base Case: Compute all the pairwise distances for each pair of points between the leaf the nodes
  void BaseCase(ANNkd_leaf* qn, ANNkd_leaf* rn) override;
  void BaseCaseCached(ANNkd_leaf* qn, ANNkd_leaf* rn) override;

  // KNN Score Function
  ANNdist Score(ANNkd_node* qn, ANNkd_node* rn) override;

  // Convert the results in the priority queue suitable for returning over to R (1-based)
  List getKnnResults();
};


#endif
