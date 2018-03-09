#include "kd_tree.h"
#include "bd_tree.h"
#include "ANN_util.h"
#include "R_kNN.h"
#include "pr_queue_k.h"
#include "node_bnb.h"

#include <Rcpp.h>
using namespace Rcpp;

#include <sys/time.h>
typedef unsigned long long timestamp_t;

static timestamp_t
  get_timestamp ()
  {
    struct timeval now;
    gettimeofday (&now, NULL);
    return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
  }

class DualTreeSearch : NodeDispatcher {	 // DTS is a node dispatcher
public:
  int n_pruned; // Performance information
  ANNkd_tree* qtree, *rtree; // query and reference tree pointers
  std::unordered_map<ANNkd_node*, node_bnb> node_info; // additional information computed at tree-construction time
  // std::set<int>* pw_check; // pairwise node-combination check
  double benchmarks[4] = { 0.0, 0.0, 0.0, 0.0};
  timestamp_t t0, t1;

  DualTreeSearch(ANNkd_tree* tree1, ANNkd_tree* tree2); // constructor

  // Pure virtuals to be overridden
  virtual void BaseCase(ANNkd_leaf* qn, ANNkd_leaf* rn) = 0;
  virtual double Score(ANNkd_node* qn, ANNkd_node* rn) = 0;
  virtual void UpdateBounds(ANNkd_leaf* n) = 0;
  virtual void UpdateBounds(ANNkd_split* n) = 0;

  // Wrapper functions
  void UpdateBounds(ANNkd_split* qn, ANNkd_split* rn) { UpdateBounds(qn); UpdateBounds(rn); }
  void UpdateBounds(ANNkd_split* qn, ANNkd_leaf* rn) { UpdateBounds(qn); UpdateBounds(rn); }
  void UpdateBounds(ANNkd_leaf* qn, ANNkd_split* rn) { UpdateBounds(qn); UpdateBounds(rn); }
  void UpdateBounds(ANNkd_leaf* qn, ANNkd_leaf* rn) { UpdateBounds(qn); UpdateBounds(rn); }

  // Default depth-first search uses dispatcher from the query tree
  void DFS(){
    // Rprintf("Starting DFS: (q:%d, r:%d)\n", getNodeID(qtree->root), getNodeID(rtree->root));
    rtree->root->ann_search_dt(1.0, *qtree->root, *this);
  }

  // For generic kd nodes
  void DFS(ANNkd_node& qn, ANNkd_node& rn){
    qn.ann_search_dt(1.0, rn, *this);
  }

  // Declare overloads for each kind of a node to dispatch
  void DFS(ANNkd_split& qn, ANNkd_split& rn);
  void DFS(ANNkd_split& qn, ANNkd_leaf& rn);
  void DFS(ANNkd_leaf& qn, ANNkd_split& rn);
  void DFS(ANNkd_leaf& qn, ANNkd_leaf& rn);

  // Various headers which may get replaced by derived implementations
  inline int getNodeID(ANNkd_node* ptr){
    return(node_info[ptr].id);
  }
  inline ANNpoint getCentroid(ANNkd_node* ptr){
    return(node_info[ptr].centroid);
  }
  inline double getMaxRadius(ANNkd_node* ptr){
    return(node_info[ptr].max_radius);
  }
  inline double getBound(ANNkd_node* ptr){
    return(node_info[ptr].bound);
  }

  // Allows setting a new bound on the node
  inline void setBound(ANNkd_node* ptr, double b){
    node_info[ptr].bound = b;
  }

  inline void printNodeBucket(ANNkd_leaf& node){
    for (int i = 0; i < node.n_pts; ++i){
      Rcout << node.bkt[i] << ", ";
    }
  }

  inline void printNodeBucket(ANNkd_split& node){}

  inline void printNodeInfo(){
    // for(auto& kv : node_info) {
    //   // List n_info = as<List>(node_info[kv.first]);
    //   // if (n_info.containsElementNamed("idx")){
    //   //   int node_id = getNodeID(kv.first);
    //   //   Rcout << "NID (" << node_id << "): ";
    //   //   printNodeBucket(*dynamic_cast<ANNkd_leaf*>(kv.first));
    //   // }
    // }
    // Rcout << "\n";
  }
};

