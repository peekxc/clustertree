#include <Rcpp.h>
using namespace Rcpp;
#include "DualTreeSearch.h"
#include "net_sort.h"

typedef unsigned int uint;
#define INDEX_TF(N,to,from) (N)*(to) - (to)*(to+1)/2 + (from) - (to) - (1)

int g = 0;

// // Initial DT search start (following triple-dispatch visitor pattern)
// void ANNkd_node::ann_search_dt(ANNdist dist, ANNkd_node& node, NodeDispatcher& dispatcher){
//   Rprintf("Called from node, recieved node\n");
//   node.ann_search_dt(dist, *this, dispatcher);
// }
// void ANNkd_node::ann_search_dt(ANNdist dist, ANNkd_split& node, NodeDispatcher& dispatcher){
//   Rprintf("Called from node, recieved split\n");
//   node.ann_search_dt(dist, *this, dispatcher);
// }
// void ANNkd_node::ann_search_dt(ANNdist dist, ANNkd_leaf& node, NodeDispatcher& dispatcher){
//   Rprintf("Called from node, recieved leaf\n");
//   node.ann_search_dt(dist, *this, dispatcher);
// }

void ANNkd_leaf::ann_search_dt(ANNdist dist, ANNkd_node& node, NodeDispatcher& dispatcher){
  // this->ann_search_dt(dist, node, dispatcher);
  // Rprintf("Called from leaf, recieved node\n");
  node.ann_search_dt(dist, *this, dispatcher);
}
void ANNkd_split::ann_search_dt(ANNdist dist, ANNkd_node& node, NodeDispatcher& dispatcher){
  //this->ann_search_dt(dist, node, dispatcher);
  // Rprintf("Called from split, recieved node\n");
  node.ann_search_dt(dist, *this, dispatcher);
}


// Dispatchers for internal nodes and leaves
void ANNkd_leaf::ann_search_dt(ANNdist dist, ANNkd_split& node, NodeDispatcher& dispatcher){
  // Rprintf("Called from leaf, recieved split\n");
  dispatcher.DFS(node, *this);
}
void ANNkd_leaf::ann_search_dt(ANNdist dist, ANNkd_leaf& node, NodeDispatcher& dispatcher){
  // Rprintf("Called from leaf, recieved leaf\n");
  dispatcher.DFS(node, *this);
}
void ANNkd_split::ann_search_dt(ANNdist dist, ANNkd_split& node, NodeDispatcher& dispatcher){
  // Rprintf("Called from split, recieved split\n");
  dispatcher.DFS(node, *this);
}
void ANNkd_split::ann_search_dt(ANNdist dist, ANNkd_leaf& node, NodeDispatcher& dispatcher){
  // Rprintf("Called from split, recieved leaf\n");
  dispatcher.DFS(node, *this);
}


// Box Decomposition Nodes (unimplemented)
void ANNbd_shrink::ann_search_dt(ANNdist dist, ANNkd_node& node, NodeDispatcher& dispatcher){}
void ANNbd_shrink::ann_search_dt(ANNdist dist, ANNkd_leaf& node, NodeDispatcher& dispatcher){}
void ANNbd_shrink::ann_search_dt(ANNdist dist, ANNkd_split& node, NodeDispatcher& dispatcher){}

// Euclidean distance
inline double computeDist(ANNpoint p1, ANNpoint p2, const int dim, ANNdist lb_threshold = ANN_DIST_INF){
  ANNdist dist = 0.0, t = 0;
  for (int i = 0; i < dim; ++i){
    t = *(p1++) - *(p2++); // compute length and adv coordinate
    // dist += (p1[i] - p2[i])*(p1[i] - p2[i]);
    // exceeds dist to k-th smallest?
    if ( (dist += ANN_POW(t)) > lb_threshold) {
      return(ANN_DIST_INF);
    }
  }
  return dist;
}

// Dual Tree Search definition
DualTreeSearch::DualTreeSearch(ANNkd_tree* tree1, ANNkd_tree* tree2) {
  qtree = tree1, rtree = tree2;
  node_info = std::unordered_map<ANNkd_node*, node_bnb>();
  g = 0;
  n_pruned = 0; // number of nodes pruned
  // base case, score, recursion, and bound update profiles
  benchmarks[0] = benchmarks[1] = benchmarks[2] = benchmarks[3] = 0.0;
}


// Declare overloads for each kind of a node to dispatch
void DualTreeSearch::DFS(ANNkd_split& qn, ANNkd_split& rn){
  // Rprintf("DFS: (q:%d [l:%d, r:%d], r:%d [l:%d, r:%d])\n",
  //         getNodeID(&qn), getNodeID(qn.child[ANN_LO]), getNodeID(qn.child[ANN_HI]),
  //         getNodeID(&rn), getNodeID(rn.child[ANN_LO]), getNodeID(rn.child[ANN_HI]));

  t0 = get_timestamp();
  double s_ll = Score(qn.child[ANN_LO], rn.child[ANN_LO]);
  double s_lh = Score(qn.child[ANN_LO], rn.child[ANN_HI]);
  double s_hl = Score(qn.child[ANN_HI], rn.child[ANN_LO]);
  double s_hh = Score(qn.child[ANN_HI], rn.child[ANN_HI]);
  t1 = get_timestamp();
  benchmarks[1] += (t1 - t0) / 1000000.0L;

  // Sort closer branches
  NodeScore scores[4] = { NodeScore(qn.child[ANN_LO], rn.child[ANN_LO], s_ll),
                          NodeScore(qn.child[ANN_LO], rn.child[ANN_HI], s_lh),
                          NodeScore(qn.child[ANN_HI], rn.child[ANN_LO], s_hl),
                          NodeScore(qn.child[ANN_HI], rn.child[ANN_HI], s_hh) } ;
  sort4_sn_stable(scores); // Sort by score using sorting network

  // Prioritize closer branches, return at first infinity. Always update the bounds after the DFS.
  // int cond = ((scores[0].score == ANN_DIST_INF ? COND_0 : 0) | (scores[1].score == ANN_DIST_INF ? COND_1 : 0) |
  //             (scores[2].score == ANN_DIST_INF ? COND_2 : 0) | (scores[3].score == ANN_DIST_INF ? COND_3 : 0));
  // switch(cond){
  //   case 0:
  //     DFS(*scores[0].lhs, *scores[0].rhs);
  //     DFS(*scores[1].lhs, *scores[1].rhs);
  //     DFS(*scores[2].lhs, *scores[2].rhs);
  //     DFS(*scores[3].lhs, *scores[3].rhs);
  //     break;
  //   case COND_3:
  //     DFS(*scores[0].lhs, *scores[0].rhs);
  //     DFS(*scores[1].lhs, *scores[1].rhs);
  //     DFS(*scores[2].lhs, *scores[2].rhs);
  //     break;
  //   case COND_2 + COND_3:
  //     DFS(*scores[0].lhs, *scores[0].rhs);
  //     DFS(*scores[1].lhs, *scores[1].rhs);
  //     break;
  //   case COND_1 + COND_2 + COND_3:
  //     DFS(*scores[0].lhs, *scores[0].rhs);
  //   case COND_0 + COND_1 + COND_2 + COND_3:
  //     return;
  //     break;
  // }
  // UpdateBounds(&qn, &rn);


  // switch(scores[0].score != ANN_DIST_INF){
  //   case true:
  //     DFS(*scores[0].lhs, *scores[0].rhs);
  //     switch(scores[1].score != ANN_DIST_INF){
  //       case true:
  //         DFS(*scores[1].lhs, *scores[1].rhs);
  //         switch(scores[2].score != ANN_DIST_INF){
  //           case true:
  //             DFS(*scores[2].lhs, *scores[2].rhs);
  //             switch(scores[3].score != ANN_DIST_INF){
  //             case true: DFS(*scores[3].lhs, *scores[3].rhs); break;
  //             default: break;
  //             }
  //             break;
  //           default:
  //             break;
  //         }
  //         break;
  //       default: break;
  //     }
  //     break;
  //   default:
  //     UpdateBounds(&qn, &rn);
  //     break;
  //
  // }
  bool sentinel = false;
  if (scores[0].score == ANN_DIST_INF) { sentinel = true; } else { DFS(*scores[0].lhs, *scores[0].rhs); }
  if (sentinel || scores[1].score == ANN_DIST_INF) { sentinel = true; } else { DFS(*scores[1].lhs, *scores[1].rhs); }
  if (sentinel || scores[2].score == ANN_DIST_INF) { sentinel = true; } else { DFS(*scores[2].lhs, *scores[2].rhs); }
  if (sentinel || scores[3].score == ANN_DIST_INF) { sentinel = true; } else { DFS(*scores[3].lhs, *scores[3].rhs); }
  UpdateBounds(&qn, &rn);
//   if (s_ll != ANN_DIST_INF){ DFS(*qn.child[ANN_LO], *rn.child[ANN_LO]); }
//   else { n_pruned++; }
//
//   if (s_lh != ANN_DIST_INF){ DFS(*qn.child[ANN_LO], *rn.child[ANN_HI]); }
//   else { n_pruned++; }
//
//   if (s_hl != ANN_DIST_INF){ DFS(*qn.child[ANN_HI], *rn.child[ANN_LO]); }
//   else { n_pruned++; }
//
//   if (s_hh != ANN_DIST_INF){ DFS(*qn.child[ANN_HI], *rn.child[ANN_HI]); }
//   else { n_pruned++; }

  // DFS(*qn.child[ANN_LO], *rn.child[ANN_LO]);
  // DFS(*qn.child[ANN_LO], *rn.child[ANN_HI]);
  // DFS(*qn.child[ANN_HI], *rn.child[ANN_LO]);
  // DFS(*qn.child[ANN_HI], *rn.child[ANN_HI]);
  // qn.child[ANN_LO]->ann_search_dt(1.0, *rn.child[ANN_LO], *this);
  // Rprintf("calling DFS: (q:%d [l:%d, r:%d], r:%d [l:%d, r:%d])\n",
  //         getNodeID(&qn), getNodeID(qn.child[ANN_LO]), getNodeID(qn.child[ANN_HI]),
  //         getNodeID(&rn), getNodeID(rn.child[ANN_LO]), getNodeID(rn.child[ANN_HI]));
  // qn.child[ANN_LO]->ann_search_dt(1.0, *rn.child[ANN_HI], *this);
  // qn.child[ANN_HI]->ann_search_dt(1.0, *rn.child[ANN_LO], *this);
  // qn.child[ANN_HI]->ann_search_dt(1.0, *rn.child[ANN_HI], *this);

};

void DualTreeSearch::DFS(ANNkd_split& qn, ANNkd_leaf& rn){
  // Rprintf("DFS: (q:%d [l:%d, r:%d], r:%d)\n",
  //         getNodeID(&qn), getNodeID(qn.child[ANN_LO]), getNodeID(qn.child[ANN_HI]),
  //         getNodeID(&rn));
  t0 = get_timestamp();
  double s_l = Score(qn.child[ANN_LO], &rn);
  double s_r = Score(qn.child[ANN_HI], &rn);
  t1 = get_timestamp();
  benchmarks[1] += (t1 - t0) / 1000000.0L;

  // Sort closer branches
  NodeScore scores[2] = { NodeScore(qn.child[ANN_LO], &rn, s_l), NodeScore(qn.child[ANN_HI], &rn, s_r) } ;
  sort2_sn_stable(scores); // Sort by score using sorting network

  // Prioritize closer branches, return at first infinity. Always update the bounds after the DFS.
  // int cond = ((scores[0].score == ANN_DIST_INF ? COND_0 : 0) | (scores[1].score == ANN_DIST_INF ? COND_1 : 0));
  // switch(cond){
  // case 0:
  //   DFS(*scores[0].lhs, *scores[0].rhs);
  //   DFS(*scores[1].lhs, *scores[1].rhs);
  //   break;
  // case COND_1:
  //   DFS(*scores[0].lhs, *scores[0].rhs);
  //   break;
  // case COND_0 + COND_1:
  //   return;
  //   break;
  // }
  // UpdateBounds(&qn, &rn);

  bool sentinel = false;
  if (scores[0].score == ANN_DIST_INF) { sentinel = true; } else { DFS(*scores[0].lhs, *scores[0].rhs); }
  if (sentinel || scores[1].score == ANN_DIST_INF) { sentinel = true; } else { DFS(*scores[1].lhs, *scores[1].rhs); }
  UpdateBounds(&qn, &rn);
 // only need to update split node
  // if (s_l != ANN_DIST_INF){ DFS(*qn.child[ANN_LO], rn); }
  // else { n_pruned++; }
  // if (s_r != ANN_DIST_INF){ DFS(*qn.child[ANN_HI], rn); }
  // else { n_pruned++; }
  // qn.child[ANN_LO]->ann_search_dt(1.0, rn, *this);
  // qn.child[ANN_HI]->ann_search_dt(1.0, rn, *this);
};

void DualTreeSearch::DFS(ANNkd_leaf& qn, ANNkd_split& rn){
  // Always call from query node
  // Rprintf("DFS: (q:%d, r:%d [l:%d, r:%d])\n",
  //         getNodeID(&qn),
  //         getNodeID(&rn), getNodeID(rn.child[ANN_LO]), getNodeID(rn.child[ANN_HI]));
  t0 = get_timestamp();
  double s_l = Score(&qn, rn.child[ANN_LO]);
  double s_r = Score(&qn, rn.child[ANN_HI]);
  t1 = get_timestamp();
  benchmarks[1] += (t1 - t0) / 1000000.0L;

  // Sort closer branches
  NodeScore scores[2] = { NodeScore(&qn, rn.child[ANN_LO], s_l), NodeScore(&qn, rn.child[ANN_HI], s_r) } ;
  sort2_sn_stable(scores); // Sort by score using sorting network

  // int cond = ((scores[0].score == ANN_DIST_INF ? COND_0 : 0) | (scores[1].score == ANN_DIST_INF ? COND_1 : 0));
  // switch(cond){
  // case 0:
  //   DFS(*scores[0].lhs, *scores[0].rhs);
  //   DFS(*scores[1].lhs, *scores[1].rhs);
  //   break;
  // case COND_1:
  //   DFS(*scores[0].lhs, *scores[0].rhs);
  //   break;
  // case COND_0 + COND_1:
  //   return;
  //   break;
  // }
  // UpdateBounds(&qn, &rn);

  bool sentinel = false;
  if (scores[0].score == ANN_DIST_INF) { sentinel = true; } else { DFS(*scores[0].lhs, *scores[0].rhs); }
  if (sentinel || scores[1].score == ANN_DIST_INF) { sentinel = true; } else { DFS(*scores[1].lhs, *scores[1].rhs); }
  UpdateBounds(&qn, &rn);
  // if (scores[0].score == ANN_DIST_INF) return; else { DFS(*scores[0].lhs, *scores[0].rhs); UpdateBounds(&rn); }
  // if (scores[1].score == ANN_DIST_INF) return; else { DFS(*scores[1].lhs, *scores[1].rhs); UpdateBounds(&rn); }
  // if (s_l != ANN_DIST_INF){ DFS(qn, *rn.child[ANN_LO]); }
  // else { n_pruned++; }
  // if (s_r != ANN_DIST_INF){ DFS(qn, *rn.child[ANN_HI]);}
  // else { n_pruned++; }
};

void DualTreeSearch::DFS(ANNkd_leaf& qn, ANNkd_leaf& rn){
  //Rcout << "(1, 1)";
  // Construct a unique index for the current pair
  // Rprintf("lDFS: (q:%d, r:%d)\n", getNodeID(&qn), getNodeID(&rn));
  // int idx = 0;
  // if (q_idx == r_idx){ idx = -(q_idx + 1); }
  // else { idx = q_idx > r_idx ? INDEX_TF(qtree->n_pts, r_idx, q_idx) : INDEX_TF(qtree->n_pts, q_idx, r_idx); }

  //printf("Checking nodes: [%d, %d] (id: %d)\n", q_idx, r_idx, idx);

  // If this specific leaf-node combination hasn't been visited before, then check base case (requires distance computations)
  // int q_idx = getNodeID(&qn), r_idx = getNodeID(&rn);
  // printf("Comparing leaf nodes: [%d, %d]\n", q_idx, r_idx);
  // if (pw_check->emplace(idx).second){ }
  BaseCase(&qn, &rn);

  // t0 = get_timestamp();
  // // UpdateBounds(&qn, &rn);
  // t1 = get_timestamp();
  // benchmarks[3] += (t1 - t0) / 1000000.0L;
}


class DualTreeKNN : public DualTreeSearch {
public:
  const int k;  // how many k-points to find the nearest neighbors of
  std::vector<ANNmin_k*> knn_pq; // knn priority queue for each point

  DualTreeKNN(ANNkd_tree* qtree, ANNkd_tree* rtree, const int k)
    : DualTreeSearch(qtree, rtree), k(k){

    // priority queue for the closest k points
    knn_pq = std::vector<ANNmin_k*>(qtree->n_pts);
    for (int i = 0; i < qtree->n_pts; ++i){ knn_pq[i] = new ANNmin_k(k); }

    // Ordered set for checking duplicate pairwise indices
    // pw_check = new set<int>();
  }

  void UpdateBounds(ANNkd_split* node){
    t0 = get_timestamp();
    double left_bound = getBound(node->child[ANN_LO]);
    double right_bound = getBound(node->child[ANN_HI]);
    double new_bound = std::max(left_bound, right_bound);
    if (new_bound > node_info[node].bound){
      node_info[node].bound = new_bound;
    }
    t1 = get_timestamp();
    benchmarks[3] += (t1 - t0) / 1000000.0L;
  }

  void UpdateBounds(ANNkd_leaf* node){
    // double max_dist = 0;
    // for (int i = 0; i < node->n_pts; ++i){
    //   int idx = node->bkt[i];
    //   double max_knn_dist = knn_pq[idx]->max_key(); // max knn dist
    //   if (max_knn_dist > max_dist){ max_dist = max_knn_dist; }
    // }
    //
    // if (max_dist > node_info[node].bound) { node_info[node].bound = max_dist; }

    // Rcout << "Bound should now be: " << getBound(node) << std::endl;
  }

  // KNN Base Case: Compute all the pairwise distances for each pair of points between the leaf the nodes
  void BaseCase(ANNkd_leaf* qn, ANNkd_leaf* rn){
    const bool same_node = qn == rn;
    t0 = get_timestamp();
    double d_pq;
    double q_max_el = 0, r_max_el = 0;
    for (int i = 0; i < qn->n_pts; ++i){
      int q_idx = qn->bkt[i];
      d_pq = knn_pq[q_idx]->max_key();
      for (int j = (same_node ? i : 0); j < rn->n_pts; ++j){ //
        int r_idx = rn->bkt[j];

        // TODO: change when qtree != rtree
        if (q_idx == r_idx) continue;
        double dist = computeDist(qtree->pts[q_idx], rtree->pts[r_idx], qtree->dim, d_pq);

        // Keep track of the maximum KNN distances
        if (dist < d_pq){
          knn_pq[q_idx]->insert(dist, r_idx);
          d_pq = knn_pq[q_idx]->max_key();
        }

        // Record up-to-date maximum knn distance
        if (d_pq > q_max_el){ q_max_el = d_pq; }

        // Also update the reference pts KNN queue, if the distance is small enough, and qtree == rtree
        // TODO: change when qtree != rtree
        if (same_node){
          double d_pr = knn_pq[r_idx]->max_key();
          if (dist < d_pr){
            knn_pq[r_idx]->insert(dist, q_idx);
            d_pr = knn_pq[r_idx]->max_key();
          }
          if (d_pr > r_max_el){ r_max_el = d_pr; }
        }
      }
    }
    t1 = get_timestamp();
    benchmarks[0] += (t1 - t0) / 1000000.0L;

    // Update max knn bound with the maximum in the priority queue
    t0 = get_timestamp();
    if (q_max_el > node_info[qn].bound) { node_info[qn].bound = q_max_el; }
    if (same_node && r_max_el > node_info[rn].bound) { node_info[rn].bound = r_max_el; }
    t1 = get_timestamp();
    benchmarks[3] += (t1 - t0) / 1000000.0L;
  }

  // KNN Score Function
  inline double Score(ANNkd_node* qn, ANNkd_node* rn){
    ANNpoint q_centroid = getCentroid(qn);
    ANNpoint r_centroid = getCentroid(rn);
    double q_bnd = getBound(qn);
    double d_min = computeDist(q_centroid, r_centroid, qtree->dim, q_bnd) - getMaxRadius(qn) - getMaxRadius(rn);
    if (d_min > 0 && d_min > q_bnd){
      // Rprintf("... pruned [%d, %d] (box_lb: %f, current bound: %f) \n", getNodeID(qn), getNodeID(rn), d_min, getBound(qn));
      return(ANN_DIST_INF); // distance between bounding boxes of nodes is high enough: prune combination
    }
    return(d_min);
  }

  // Convert the results in the priority queue suitable for returning over to R (1-based)
  List getKnnResults(){
    const int n = qtree->n_pts;
    IntegerMatrix idx = Rcpp::no_init_matrix(n, k);
    NumericMatrix dist = Rcpp::no_init_matrix(n, k);
    for (int i = 0; i < n; ++i){
      for (int k_i = 0; k_i < k; ++k_i){
        dist(i, k_i) = (ANNdist) knn_pq[i]->ith_smallest_key(k_i);
        idx(i, k_i) = (int) knn_pq[i]->ith_smallest_info(k_i);
      }
    }
    return(List::create(_["dist"] = dist, _["idx"] = idx + 1));
  }
};




// Iterate in a depth-first search manner through the tree to collect all of the nodes
// void DFS(ANNkd_leaf* cnode, std::vector<ANNkd_node*>& all_nodes){
//   all_nodes.push_back(cnode);
// }
// void DFS(ANNkd_split* cnode, std::vector<ANNkd_node*>& all_nodes){
//   all_nodes.push_back(cnode);
//   DFS(cnode->child[ANN_LO], all_nodes);
//   DFS(cnode->child[ANN_HI], all_nodes);
// }

// [[Rcpp::export]]
void testTrees(const NumericMatrix& x, const int k, const int bucketSize = 15){
  // Copy data over to ANN point array
  ANNpointArray queryPts = matrixToANNpointArray(x);

  // Test creating a regular tree
  timestamp_t t0 = get_timestamp();
  ANNkd_tree* kdTree1 = new ANNkd_tree(queryPts, x.nrow(), x.ncol(), bucketSize, (ANNsplitRule) 5);
  timestamp_t t1 = get_timestamp();
  Rcout << "Seconds to create regular tree: " << (t1 - t0) / 1000000.0L << std::endl;

  // Test creating a dual-tree
  t0 = get_timestamp();
  std::unordered_map<ANNkd_node*, node_bnb> node_info; // container to store various properties related to tree construction
  ANNkd_tree* kdTree2 = new ANNkd_tree(queryPts, x.nrow(), x.ncol(), bucketSize, (ANNsplitRule) 5, node_info);
  t1 = get_timestamp();
  Rcout << "Seconds to create dual tree: " << (t1 - t0) / 1000000.0L << std::endl;

  delete kdTree1;
  delete kdTree2;
}


// [[Rcpp::export]]
List testDTS(const NumericMatrix& x, const int k, const int bucketSize = 15){

  timestamp_t t0, t1;

  // Copy data over to ANN point array
  t0 = get_timestamp();
  ANNpointArray queryPts = matrixToANNpointArray(x);
  t1 = get_timestamp();
  Rcout << "Seconds to copy points: " << (t1 - t0) / 1000000.0L << std::endl;

  // Create the dual tree
  t0 = get_timestamp();
  std::unordered_map<ANNkd_node*, node_bnb> node_info; // container to store various properties related to tree construction
  ANNkd_tree* kdTreeR = new ANNkd_tree(queryPts, x.nrow(), x.ncol(), bucketSize, (ANNsplitRule) 5, node_info);
  t1 = get_timestamp();
  Rcout << "Seconds to create tree: " << (t1 - t0) / 1000000.0L << std::endl;

  // Create the Dual Tree KNN object, pointing to the new node info
  t0 = get_timestamp();
  DualTreeKNN dts = DualTreeKNN(kdTreeR, kdTreeR, k);
  dts.node_info = node_info;
  t1 = get_timestamp();
  Rcout << "Seconds to create DualTree object: " << (t1 - t0) / 1000000.0L << std::endl;

  // for(auto& kv : dts.node_info) {
  //   Rcout << dts.getBound(kv.first) << std::endl;
  // }

  // Depth-First Traversal solves the KNN
  timestamp_t t0_cpy = get_timestamp();
  dts.DFS();
  timestamp_t t1_cpy = get_timestamp();
  double total_time = (t1_cpy - t0_cpy) / 1000000.0L;
  Rcout << "Seconds to do DFS: " << total_time << std::endl;

  // base case, score, recursion, and bound update profiles
  Rprintf("Benchmarks: [ Base case: %f, Score: %f, Recursion: %f, Bound Update: %f, Unaccounted for: %f]\n",
          dts.benchmarks[0], dts.benchmarks[1], dts.benchmarks[2], dts.benchmarks[3],
          total_time-(dts.benchmarks[0]+dts.benchmarks[1]+dts.benchmarks[2]+dts.benchmarks[3]));

  // dts.printNodeInfo();
  /// Rcout << "Node combinations pruned: " << dts.n_pruned << std::endl;

  // ANNkdStats tree_stat = ANNkdStats();
  // kdTreeR->getStats(tree_stat);
  // Rcout << "Tree Depth: " << tree_stat.depth << std::endl;
  // Rcout << "Number of Split Nodes: " << tree_stat.n_spl << std::endl;

  // Retrieve all properties built from the tree as a List
  // List tree_properties = List();
  // for(auto& kv : dts.node_info) {
  //   tree_properties.push_back(kv.second);
  //   double mr = dts.getMaxRadius(kv.first);
  //   ANNpoint pt = dts.getCentroid(kv.first);
  //   //Rcout << mr << std::endl;
  //   //return(tree_properties);
  // }
  // Rcout << tree_properties.size() << std::endl;
  t0 = get_timestamp();
  List res = List::create(_["knn"] = dts.getKnnResults());
  t1 = get_timestamp();
  Rcout << "Seconds to copy results: " << (t1 - t0) / 1000000.0L << std::endl;
  return(res);
  // return(List::create(_["knn"] = dts.getKnnResults(),
  //                     _["info"] = tree_properties));
}

/*** R
# x <- as.matrix(iris[, 1:4])
set.seed(1234)
x <- rbind(cbind(rnorm(5), rnorm(5)),
           cbind(rnorm(5, mean = 5), rnorm(5, mean = 5)))


info <- clustertree:::testDTS(x, k = 5, bucketSize = 1) ## 11,19 should be pruned
info_db <- dbscan::kNN(x, k = 5)
all((info_db$id - 1) == info$knn$idx)
all((info_db$dist - sqrt(info$knn$dist)) < .Machine$double.eps)

n <- 5000
x2 <- rbind(cbind(rnorm(n), rnorm(n), rnorm(n)),
            cbind(rnorm(n, mean = 5), rnorm(n, mean = 5), rnorm(n, mean = 5)))
all(dim(unique(x2)) == dim(x2))
## Benchmarks: [ Base case: 0.906046, Score: 0.682107, Recursion: 0.000000, Bound Update: 0.635919 ] (n == 5000, k == 25, bucketSize == 10)
info <- clustertree:::testDTS(x2, k = 25, bucketSize = 10) ## 11,19 should be pruned
info_db <- dbscan::kNN(x2, k = 25, bucketSize = 10)
all(info_db$id == info$knn$idx)
all((info_db$dist - sqrt(info$knn$dist)) < .Machine$double.eps)

## 2.27183
## Testing tree construction time
clustertree:::testTrees(x2, k = 25, bucketSize = 10) ## 11,19 should be pruned

## Profiling
microbenchmark::microbenchmark({ dbscan::kNN(x2, k = 25, bucketSize = 10) }, times = 1L)
microbenchmark::microbenchmark({ clustertree:::testDTS(x2, k = 25, bucketSize = 10) }, times = 1L)

invisible(gprofiler::profile({ replicate(5, invisible(clustertree:::testDTS(x2, k = 25, bucketSize = 10))) }, filename = "cluster_perf.out"))
system(command = "pprof.pl --web /Library/Frameworks/R.framework/Resources/bin/exec/R /Users/mpiekenbrock/clustertree/clustertree_prof.out")

plotLeafBB <- function(nodes, leaf_only = FALSE){
  for (node in nodes){
    lines(x = c(node$lb[1], node$ub[1]), y = c(node$lb[2], node$lb[2]))
    lines(x = c(node$ub[1], node$ub[1]), y = c(node$lb[2], node$ub[2]))
    lines(x = c(node$ub[1], node$lb[1]), y = c(node$ub[2], node$ub[2]))
    lines(x = c(node$lb[1], node$lb[1]), y = c(node$ub[2], node$lb[2]))
    points(node$centroid[[1]], node$centroid[[2]], col = "red", pch = 20, cex = 0.5)
    plotrix::draw.circle(x = node$centroid[[1]], y = node$centroid[[2]], radius = node$max_radius)
  }
}
plot(x, xlim = range(x[, 1]) + c(-1, 1), ylim = range(x[, 2]) + c(-1, 1), asp = 1)
text(x, labels = 0:(nrow(x) - 1), pos = 3)
node_idx <- sapply(info$info, function(node) node$id)
plotLeafBB(info$info[which(node_idx %in% c(11, 12))])


which(sapply(sapply(info$info, function(node) node$id), function(id) if(is.null(id)) FALSE else id == 1L))
score(info$info[[8]], info$info[[2]])


score <- function(qn, rn){
  dist_qr <- dist(rbind(qn$centroid, rn$centroid))[[1]]^2
  dist_qr - qn$max_radius - rn$max_radius
}

xb <- cbind(rnorm(1000), rnorm(1000)) #as.matrix(iris[, 1:4]) #
cl1 <- function(){ invisible(clustertree:::testDTS(xb, k = 15, bucketSize = 15L)) }
cl2 <- function(){ invisible(dbscan::kNN(xb, k = 15, bucketSize = 15L)) }
microbenchmark::microbenchmark(cl1(), times = 3L)
microbenchmark::microbenchmark(cl2(), times = 3L)

all((cl2()$idx - 1) == cl1()$idx)

*/
