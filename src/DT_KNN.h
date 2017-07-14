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

// Use a fixed-sized scoremap to easily score and keep track of which recursions to prioritize
template <typename T, typename V> struct KeyValue {
  T key;
  V value;
  KeyValue(){};
  KeyValue(T _key, V _val) : key(_key), value(_val){}
  bool operator<(KeyValue kv) const
  { return value > kv.value; }
};
struct ScoreMap {
  std::array< KeyValue<ANNdist, int>, 4> scores;
  std::array< NODE_PAIR, 4 > nodes;
  ScoreMap() {
    scores = std::array< KeyValue<ANNdist, int>, 4>();
    scores[0].key = 0, scores[0].value = ANN_DIST_INF;
    scores[1].key = 1, scores[1].value = ANN_DIST_INF;
    scores[2].key = 2, scores[2].value = ANN_DIST_INF;
    scores[3].key = 3, scores[3].value = ANN_DIST_INF;
    nodes = std::array< NODE_PAIR, 4 >();
  }
  void inline reset(){
    scores[0].key = 0, scores[0].value = ANN_DIST_INF;
    scores[1].key = 1, scores[1].value = ANN_DIST_INF;
    scores[2].key = 2, scores[2].value = ANN_DIST_INF;
    scores[3].key = 3, scores[3].value = ANN_DIST_INF;
  }
  void inline insert(NODE_PAIR node_pair, ANNdist score, int idx){
    nodes[idx] = node_pair;
    scores[idx].key = score, scores[idx].value = idx;
  }
  inline const NODE_PAIR& operator[](size_t idx) const { return nodes[idx]; }
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
    R_INFO("Min dist. Q (" << N_q << ") <--> R (" << N_q << "): " << min_dist_qr << " (smaller than? " << max_knn_dist << " )\n")

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

    // Compute distance and retrieve KNN pr queue
    ANNdist dist = annDist(d, qtree->pts[q_idx], rtree->pts[r_idx]);// (ANNdist) ANN_SUM(box_dist, ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));
    ANNmin_k& q_knn = *knn->at(q_idx); // Retrieve query point's knn priority queue

    // Is the computed distance less than the current k-nearest neighbor?
    if (dist < q_knn.max_key()){
      q_knn.insert(dist, r_idx);
      // R_INFO("Storing better neighbor for: " << q_idx << " (id=" << r_idx << ", dist=" << dist << ")\n")
    }

    // If the query == reference sets, set the other nearest neighbor as well!
    if (knn_identity && q_idx != r_idx && dist < knn->at(r_idx)->max_key()){
      knn->at(r_idx)->insert(dist, q_idx);
      // R_INFO("Storing better neighbor for: " << r_idx << " (id=" << q_idx << ", dist=" << dist << ")\n")
    }
  };
};