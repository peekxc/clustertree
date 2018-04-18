#ifndef DUALTREESEARCH_H
#define DUALTREESEARCH_H

#include "RcppHeader.h"
#include "kd_tree.h"
#include "bd_tree.h"
#include "ANN_util.h"
#include "R_kNN.h"
#include "pr_queue_k.h"
#include "node_bnb.h"
#include "kd_util.h"
#include "net_sort.h"
#include "simple_structs.h"
#include "profiler.h"
#include "metrics.h"

// Cut differences enumeration. CUT_LO indicates that a query nodes convex subset along dimension d lies entirely
// lower / 'to the left' of the reference nodes. CUT_HI is defined inversely. CUT_OVERLAP is reserved for the case where
// the intervals between the query and reference nodes overlap in d.
enum CUT_DIFF { CUT_LO = -1, CUT_OVERLAP = 0, CUT_HI = 1 };

// Quicker hash function to improve unordered_map access
template <class T>
inline void hash_combine(std::size_t & seed, const T & v)
{
  std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}
template<typename S, typename T> struct hash<pair<S, T>>
{
  inline size_t operator()(const pair<S, T> & v) const
  {
    size_t seed = 0;
    ::hash_combine(seed, v.first);
    ::hash_combine(seed, v.second);
    return seed;
  }
};

struct cache {
  bool reverse_ordered;
  std::vector<ANNdist>* dist_map;
  cache() : reverse_ordered(false), dist_map(nullptr){ }
  cache(bool order, std::vector<ANNdist>* dist_cache) : reverse_ordered(order), dist_map(dist_cache){ }
};
typedef std::unordered_map< pair<ANNkd_node*, ANNkd_node*>, cache > CACHE;


template <typename METRIC_T>
class DualTreeSearch : NodeDispatcher {	 // DTS is a node dispatcher
public:
  static_assert(std::is_base_of<Metric, METRIC_T>::value, "Must be a valid metric.");
  METRIC_T metric; // metric distance to use, includes dimensionality of space
  ANNkd_tree* qtree, *rtree; // query and reference tree pointers
  NODE_INFO qinfo, rinfo; // additional node information computed at tree-construction time

  std::vector<CUT_DIFF> cut_dir; // cut directions; one of -1
  // Performance information
  int n_pruned; // number of branches pruned
  int n_comparisons; // number of point comparisons
  int n_accesses; // number of accesses to the dist map

  // Auxilliary containers
  CACHE cache_map; // optional cache for already computed distances

  // Constructor
  DualTreeSearch(ANNkd_tree* qtree, ANNkd_tree* rtree, const NODE_INFO& qinfo, const NODE_INFO& rinfo, METRIC_T& m);

  // Pure virtuals to be overridden by derived classes
  virtual void BaseCase(ANNkd_leaf* qn, ANNkd_leaf* rn) = 0;
  virtual void BaseCaseCached(ANNkd_leaf* qn, ANNkd_leaf* rn) = 0; // cached version
  virtual ANNdist Score(ANNkd_node* qn, ANNkd_node* rn) = 0;
  virtual void UpdateBounds(ANNkd_leaf* n, NODE_INFO& ninfo) = 0;
  virtual void UpdateBounds(ANNkd_split* n, NODE_INFO& ninfo) = 0;

  // Wrapper functions
  void UpdateBounds(ANNkd_split* qn, ANNkd_split* rn) { UpdateBounds(qn, qinfo); UpdateBounds(rn, rinfo); }
  void UpdateBounds(ANNkd_split* qn, ANNkd_leaf* rn) { UpdateBounds(qn, qinfo); UpdateBounds(rn, rinfo); }
  void UpdateBounds(ANNkd_leaf* qn, ANNkd_split* rn) { UpdateBounds(qn, qinfo); UpdateBounds(rn, rinfo); }
  void UpdateBounds(ANNkd_leaf* qn, ANNkd_leaf* rn) { UpdateBounds(qn, qinfo); UpdateBounds(rn, rinfo); }
  void resetBounds();

  // Default depth-first search uses dispatcher from the query tree
  void DFS();

  // Declare overloads for each kind of a node to dispatch
  void DFS(ANNkd_node& qn, ANNkd_node& rn) { qn.ann_search_dt(rn, *this); };
  void DFS(ANNkd_split& qn, ANNkd_split& rn) override;
  void DFS(ANNkd_split& qn, ANNkd_leaf& rn) override;
  void DFS(ANNkd_leaf& qn, ANNkd_split& rn) override;
  void DFS(ANNkd_leaf& qn, ANNkd_leaf& rn) override;

  /* === Various inlines for accessing node information === */

  // Query node information
  inline int getQN_ID(ANNkd_node* ptr){
    if (ptr->id > qinfo.size()){
      R_OUT("Q id: %d greater than id size %d\n", ptr->id, qinfo.size());
      stop("here");
    }
    return(qinfo[ptr->id].id);
  }
  inline ANNpoint getQN_Centroid(ANNkd_node* ptr){ return(qinfo[ptr->id].centroid); }
  inline double getQN_MaxRadius(ANNkd_node* ptr){ return(qinfo[ptr->id].max_radius); }
  inline double getQN_Bound(ANNkd_node* ptr){ return(qinfo[ptr->id].bound); }
  inline void setQN_Bound(ANNkd_node* ptr, double b){ qinfo[ptr->id].bound = b; }

  // Reference node information
  inline int getRN_ID(ANNkd_node* ptr){
    if (ptr->id > rinfo.size()){
      R_OUT("R id: %d greater than id size %d\n", ptr->id, rinfo.size());
      stop("here");
    }
    return(rinfo[ptr->id].id);
  }
  inline ANNpoint getRN_Centroid(ANNkd_node* ptr){ return(rinfo[ptr->id].centroid); }
  inline double getRN_MaxRadius(ANNkd_node* ptr){ return(rinfo[ptr->id].max_radius); }
  inline double getRN_Bound(ANNkd_node* ptr){ return(rinfo[ptr->id].bound); }
  inline void setRN_Bound(ANNkd_node* ptr, double b){ rinfo[ptr->id].bound = b; }

  // // Generic function to compute distance
  // inline double computeDist(ANNpoint p1, ANNpoint p2, const int dim, ANNdist lb_threshold = ANN_DIST_INF) {
  //   ANNdist dist = 0.0, t = 0;
  //   for (int i = 0; i < dim; ++i){
  //     t = *(p1++) - *(p2++); // compute length and adv coordinate
  //     // exceeds dist to k-th smallest?
  //     if ( (dist = ANN_SUM(dist, ANN_POW(t))) > lb_threshold) {
  //       return(ANN_DIST_INF);
  //     }
  //   }
  //   return dist;
  // }

  // Computes or retrieves a precomputed pairwise distance between points
  ANNdist getPairDist(const int q_idx, const int r_idx, bool use_cache = true);


  // --- Bounding functions ---
  // (Optimal) box distance can be used
  // template<class M, class std::enable_if<!std::is_base_of<MinkowskiMetric, M>::value, ANNdist>::type = 0>
  double computeBound(ANNkd_node* qn, ANNkd_node* rn);

  // Circumsphere bound, if box distance not available
  // template<class M, class std::enable_if<std::is_base_of<MinkowskiMetric, M>::value, ANNdist>::type = 0>
  // double computeBound(ANNkd_node* qn, ANNkd_node* rn);

  void updateCutDirections(ANNkd_node* qn, ANNkd_node* rn);
  ANNdist computeBoxDist(ANNkd_node* qn, ANNkd_node* rn);

  // Utility
  List getTreeProperties(bool query = true);
  inline void printNodeBucket(ANNkd_leaf& node){
    for (int i = 0; i < node.n_pts; ++i){
      Rcout << node.bkt[i] << ", ";
    }
  }
  inline void printNodeBucket(ANNkd_split& node){}
  inline void printNodeInfo(){}
};

#endif
