#include "DualTreeSearch.h"

typedef unsigned int uint;
#define INDEX_TF(N,to,from) (N)*(to) - (to)*(to+1)/2 + (from) - (to) - (1)

#define ANN_PERF

// Initial DT search start (following triple-dispatch visitor pattern)
void ANNkd_leaf::ann_search_dt(ANNkd_node& node, NodeDispatcher& dispatcher){
  // this->ann_search_dt(dist, node, dispatcher);
  // Rprintf("Called from leaf, recieved node\n");
  node.ann_search_dt(*this, dispatcher);
}
void ANNkd_split::ann_search_dt(ANNkd_node& node, NodeDispatcher& dispatcher){
  //this->ann_search_dt(dist, node, dispatcher);
  // Rprintf("Called from split, recieved node\n");
  node.ann_search_dt(*this, dispatcher);
}


// Dispatchers for internal nodes and leaves
void ANNkd_leaf::ann_search_dt(ANNkd_split& node, NodeDispatcher& dispatcher){
  // Rprintf("Called from leaf, recieved split\n");
  dispatcher.DFS(node, *this);
}
void ANNkd_leaf::ann_search_dt(ANNkd_leaf& node, NodeDispatcher& dispatcher){
  // Rprintf("Called from leaf, recieved leaf\n");
  dispatcher.DFS(node, *this);
}
void ANNkd_split::ann_search_dt(ANNkd_split& node, NodeDispatcher& dispatcher){
  // Rprintf("Called from split, recieved split\n");
  dispatcher.DFS(node, *this);
}
void ANNkd_split::ann_search_dt(ANNkd_leaf& node, NodeDispatcher& dispatcher){
  // Rprintf("Called from split, recieved leaf\n");
  dispatcher.DFS(node, *this);
}

// Box Decomposition Nodes (unimplemented)
void ANNbd_shrink::ann_search_dt(ANNkd_node& node, NodeDispatcher& dispatcher){}
void ANNbd_shrink::ann_search_dt(ANNkd_leaf& node, NodeDispatcher& dispatcher){}
void ANNbd_shrink::ann_search_dt(ANNkd_split& node, NodeDispatcher& dispatcher){}// Euclidean distance


// Reset the bounds associated with each node
template <class METRIC_T>
void DualTreeSearch<METRIC_T>::resetBounds(){
  std::for_each(qinfo.begin(), qinfo.end(), [&](node_bnb& ni){ ni.bound = ANN_DIST_INF; });
  if (std::addressof(qinfo) != std::addressof(rinfo)){ std::for_each(rinfo.begin(), rinfo.end(), [&](node_bnb& ni){ ni.bound = ANN_DIST_INF; }); }
};

// Dual Tree Search definition
template <class METRIC_T>
DualTreeSearch<METRIC_T>::DualTreeSearch(ANNkd_tree* queryTree, ANNkd_tree* refTree, const NODE_INFO& queryInfo, const NODE_INFO& refInfo, METRIC_T& m) : metric(m) {
  qtree = queryTree, rtree = refTree;
  qinfo = queryInfo, rinfo = refInfo;
  cut_dir = std::vector<CUT_DIFF>(qtree->dim);  // Vector to record the cut differences
  cache_map = CACHE();

  // if logging is enabled
  FILE * pFile = fopen ("log.txt","w");
  if (pFile!=NULL) { fclose (pFile); }
}

// Default to start the DFS
template <class METRIC_T>
void DualTreeSearch<METRIC_T>::DFS(){
  // Start the depth first search with the root nodes
  R_OUT("Starting DFS: (%d, %d)\n", qinfo[qtree->root->id].id, rinfo[rtree->root->id].id);
  // rtree->root->ann_search_dt(*qtree->root, *this);
  qtree->root->ann_search_dt(*rtree->root, *this);
}

// Function to update the cut directions
// The cut direction maintains which coordinates of the query nodes convex subset
// lie closest to the reference nodes convex subset, i.e. either whether the lower/upper
// bound is closer (-1 and 1), or if the intervals intersect (0).
template <class METRIC_T>
void DualTreeSearch<METRIC_T>::updateCutDirections(ANNkd_node* qn, ANNkd_node* rn){
  // for (int d_i = 0; d_i < qtree->dim; ++d_i){
  //   if (node_info[qn->id].hi[d_i] < node_info[rn->id].lo[d_i]){
  //     cut_dir[d_i] = CUT_LO;
  //   } else if (node_info[qn->id].lo[d_i] > node_info[rn->id].hi[d_i]){
  //     cut_dir[d_i] = CUT_HI;
  //   } else {
  //     cut_dir[d_i] = CUT_OVERLAP;
  //   }
  // }
}

// Bounding functions; If the metric is a minkowski metric, the (optimal) box distance can be used, otherwise the
// circumsphere bound may be used if the coordinates cannot be embedded in the Cartesian plane.
template<class METRIC_T>
double DualTreeSearch<METRIC_T>::computeBound(ANNkd_node* qn, ANNkd_node* rn){

    // Using the actual box distance as the bound on the smallest kNN distance
    ADD_BENCHMARK_DESCRIPTION("Updating cut directions", 11)
    START_TIME()
    // updateCutDirections(qn, rn); // Update which coordinates will need to compute box distance
    END_TIME(11)

    ADD_BENCHMARK_DESCRIPTION("Computing box distance", 12)
    START_TIME()
    ANNdist actual_box_dist = computeBoxDist(qn, rn); // relies on updated cut directions!
    END_TIME(12)
    return(actual_box_dist);

    // Using the circumsphere bound
    ANNpoint q_centroid = qinfo[qn->id].centroid;
    ANNpoint r_centroid = rinfo[rn->id].centroid;
    ANNdist bound = metric.finalize(metric.distance(q_centroid, r_centroid)) - (qinfo[qn->id].max_radius + rinfo[rn->id].max_radius);
    // R_OUT("Centroid dist: %f, max qn: %f, max rn: %f\n", annDist(qtree->dim, q_centroid, r_centroid), std::pow(getMaxRadius(qn), 2), std::pow(getMaxRadius(rn), 2));
    return(bound < 0 ? 0 : std::pow(bound, 2));

  }

// Using the circumsphere bound
// template<class METRIC_T>
// double DualTreeSearch<METRIC_T>::computeBound(ANNkd_node* qn, ANNkd_node* rn){
//
// }


// Declare overloads for each kind of a node to dispatch
template <class METRIC_T>
void DualTreeSearch<METRIC_T>::DFS(ANNkd_split& qn, ANNkd_split& rn){
  { R_OUT("actual DFS: (q:%d [l:%d, r:%d], r:%d [l:%d, r:%d])\n",
          qn.id, qn.child[ANN_LO]->id, qn.child[ANN_HI]->id,
          rn.id, rn.child[ANN_LO]->id, rn.child[ANN_HI]->id); }
  R_OUT("DFS: (q:%d [l:%d, r:%d], r:%d [l:%d, r:%d])\n",
        getQN_ID(&qn), getQN_ID(qn.child[ANN_LO]), getQN_ID(qn.child[ANN_HI]),
        getRN_ID(&rn), getRN_ID(rn.child[ANN_LO]), getRN_ID(rn.child[ANN_HI]));

  // ADD_BENCHMARK_DESCRIPTION("Updating cut directions", 11)
  // START_TIME()
  // const int q_d_i = qn.cut_dim, r_d_i = rn.cut_dim;
  // cut_dir[q_d_i] = node_info[&qn].hi[q_d_i] < node_info[&rn].lo[q_d_i] ? CUT_LO : node_info[&qn].lo[q_d_i] > node_info[&rn].hi[q_d_i] ? CUT_HI : CUT_OVERLAP;
  // cut_dir[r_d_i] = node_info[&qn].hi[r_d_i] < node_info[&rn].lo[r_d_i] ? CUT_LO : node_info[&qn].lo[r_d_i] > node_info[&rn].hi[r_d_i] ? CUT_HI : CUT_OVERLAP;
  // END_TIME(11)

  ADD_BENCHMARK_DESCRIPTION("Scoring nodes, total", 19)
  START_NEW_TIME()
  double s_ll = Score(qn.child[ANN_LO], rn.child[ANN_LO]);
  double s_lh = Score(qn.child[ANN_LO], rn.child[ANN_HI]);
  double s_hl = Score(qn.child[ANN_HI], rn.child[ANN_LO]);
  double s_hh = Score(qn.child[ANN_HI], rn.child[ANN_HI]);
  END_NEW_TIME(19)


  // Sort closer branches
  ADD_BENCHMARK_DESCRIPTION("Storing and sorting node scores", 17)
  START_TIME()
  NodeScore scores[4] = { NodeScore(qn.child[ANN_LO], rn.child[ANN_LO], s_ll),
                          NodeScore(qn.child[ANN_LO], rn.child[ANN_HI], s_lh),
                          NodeScore(qn.child[ANN_HI], rn.child[ANN_LO], s_hl),
                          NodeScore(qn.child[ANN_HI], rn.child[ANN_HI], s_hh) } ;
  sort4_sn_stable(scores); // Sort by score using sorting network
  END_TIME(17)

  const int n_branches = 4;
  for (int i = 0; i < n_branches; ++i){
    if (scores[i].score != ANN_DIST_INF){
      DFS(*scores[i].lhs, *scores[i].rhs);
      ADD_BENCHMARK_DESCRIPTION("Updating split node bounds", 13)
      START_TIME()
      R_OUT("Updating bounds--> q:%d, r:%d\n", qn.id, rn.id);
      UpdateBounds(&qn, &rn);
      END_TIME(13)
    } else {
      break; // break out of loop and end recursion at first infinity
    }
  }
  return;
};

// Testing triple recursion instead of quad:
// void DualTreeSearch::DFS(ANNkd_split& qn, ANNkd_split& rn){
//   R_OUT("DFS: (q:%d [l:%d, r:%d], r:%d [l:%d, r:%d])\n",
//           getNodeID(&qn), getNodeID(qn.child[ANN_LO]), getNodeID(qn.child[ANN_HI]),
//           getNodeID(&rn), getNodeID(rn.child[ANN_LO]), getNodeID(rn.child[ANN_HI]));
//
//   double s_ll = Score(qn.child[ANN_LO], rn.child[ANN_LO]);
//   double s_lh = Score(qn.child[ANN_LO], rn.child[ANN_HI]);
//   double s_hh = Score(qn.child[ANN_HI], rn.child[ANN_HI]);
//   if (&qn != &rn){
//     double s_hl = Score(qn.child[ANN_HI], rn.child[ANN_LO]);
//     // Sort closer branches
//     NodeScore scores[4] = { NodeScore(qn.child[ANN_LO], rn.child[ANN_LO], s_ll, LL),
//                             NodeScore(qn.child[ANN_LO], rn.child[ANN_HI], s_lh, LR),
//                             NodeScore(qn.child[ANN_HI], rn.child[ANN_LO], s_lh, RL),
//                             NodeScore(qn.child[ANN_HI], rn.child[ANN_HI], s_hh, RR) } ;
//     sort4_sn_stable(scores); // Sort by score using sorting network
//     const int n_branches = 4;
//     for (int i = 0; i < n_branches; ++i){
//       if (scores[i].score != ANN_DIST_INF){
//         DFS(*scores[i].lhs, *scores[i].rhs);
//         UpdateBounds(&qn, &rn);
//       } else {
//         break; // break out of loop and end recursion at first infinity
//       }
//     }
//   } else {
//     // Sort closer branches
//     NodeScore scores[3] = { NodeScore(qn.child[ANN_LO], rn.child[ANN_LO], s_ll, LL),
//                             NodeScore(qn.child[ANN_LO], rn.child[ANN_HI], s_lh, LR),
//                             NodeScore(qn.child[ANN_HI], rn.child[ANN_HI], s_hh, RR) } ;
//     sort3_sn_stable(scores); // Sort by score using sorting network
//     const int n_branches = 3;
//     for (int i = 0; i < n_branches; ++i){
//       if (scores[i].score != ANN_DIST_INF){
//         DFS(*scores[i].lhs, *scores[i].rhs);
//         UpdateBounds(&qn, &rn);
//       } else {
//         break; // break out of loop and end recursion at first infinity
//       }
//     }
//   }
//   return;
// };

template <class METRIC_T>
void DualTreeSearch<METRIC_T>::DFS(ANNkd_split& qn, ANNkd_leaf& rn){
  { R_OUT("actual DFS: (q:%d [l:%d, r:%d], r:%d)\n", qn.id, qn.child[ANN_LO]->id, qn.child[ANN_HI]->id, rn.id); }
  R_OUT("DFS: (q:%d [l:%d, r:%d], r:%d)\n",
        getQN_ID(&qn), getQN_ID(qn.child[ANN_LO]), getQN_ID(qn.child[ANN_HI]),
        getRN_ID(&rn));

  // ADD_BENCHMARK_DESCRIPTION("Updating cut directions", 11)
  // START_TIME()
  // const int d_i = qn.cut_dim;
  // cut_dir[d_i] = node_info[&qn].hi[d_i] < node_info[&rn].lo[d_i] ? CUT_LO : node_info[&qn].lo[d_i] > node_info[&rn].hi[d_i] ? CUT_HI : CUT_OVERLAP;
  // END_TIME(11)

  ADD_BENCHMARK_DESCRIPTION("Scoring nodes, total", 19)
  START_NEW_TIME()
  double s_l = Score(qn.child[ANN_LO], &rn);
  double s_r = Score(qn.child[ANN_HI], &rn);
  END_NEW_TIME(19)

  // Sort closer branches
  ADD_BENCHMARK_DESCRIPTION("Storing and sorting node scores", 17)
  START_TIME()
  NodeScore scores[2] = { NodeScore(qn.child[ANN_LO], &rn, s_l), NodeScore(qn.child[ANN_HI], &rn, s_r) } ;
  sort2_sn_stable(scores); // Sort by score using sorting network
  END_TIME(17)

  // Prioritize closer branches, return at first infinity. Always update the bounds after the DFS.
  const int n_branches = 2;
  for (int i = 0; i < n_branches; ++i){
    if (scores[i].score != ANN_DIST_INF){
      DFS(*scores[i].lhs, *scores[i].rhs);
      ADD_BENCHMARK_DESCRIPTION("Updating split node bounds", 13)
      START_TIME()
      R_OUT("Updating bounds--> q:%d, r:%d\n", qn.id, rn.id);
      UpdateBounds(&qn, &rn);
      END_TIME(13)
    } else {
      break; // break out of loop and end recursion at first infinity
    }
  }
  return;
  // int cond = ((scores[0].score == ANN_DIST_INF ? COND_0 : 0) | (scores[1].score == ANN_DIST_INF ? COND_1 : 0));
  // switch(cond){
  // case 0:
  //   R_OUT("Base cond: Recurse all\n", 0);
  //   PathDFS(lpath, rpath, scores[0]), UpdateBounds(&qn);
  //   // R_OUT("-- Continuing from parent: (q:%d, r:%d) --> q_: %d, r_: %d\n", getNodeID(&qn), getNodeID(&rn), getNodeID(scores[1].lhs), getNodeID(scores[1].rhs));
  //   PathDFS(lpath, rpath, scores[1]), UpdateBounds(&qn);
  //   break;
  // case COND_1:
  //   R_OUT("Cond 1: Recurse 0 (%d, %d), Prune 1 (%d, %d)\n", scores[0].lhs, scores[0].rhs, scores[1].lhs, scores[1].rhs);
  //   PathDFS(lpath, rpath, scores[0]), UpdateBounds(&qn); break;
  // case COND_0 + COND_1: return; break;
  // }
  // // UpdateBounds(&qn); // Only need to update split node
};

template <class METRIC_T>
void DualTreeSearch<METRIC_T>::DFS(ANNkd_leaf& qn, ANNkd_split& rn){
  { R_OUT("actual DFS: (q:%d, r:%d [l:%d, r:%d])\n", qn.id, rn.id, rn.child[ANN_LO]->id, rn.child[ANN_HI]->id); }
  R_OUT("DFS: (q:%d, r:%d [l:%d, r:%d])\n",
        getQN_ID(&qn),
        getRN_ID(&rn), getRN_ID(rn.child[ANN_LO]), getRN_ID(rn.child[ANN_HI]));

  // ADD_BENCHMARK_DESCRIPTION("Updating cut directions", 11)
  // START_TIME()
  // const int d_i = rn.cut_dim;
  // cut_dir[d_i] = node_info[&qn].hi[d_i] < node_info[&rn].lo[d_i] ? CUT_LO : node_info[&qn].lo[d_i] > node_info[&rn].hi[d_i] ? CUT_HI : CUT_OVERLAP;
  // END_TIME(11)

  ADD_BENCHMARK_DESCRIPTION("Scoring nodes, total", 19)
  START_NEW_TIME()
  double s_l = Score(&qn, rn.child[ANN_LO]);
  double s_r = Score(&qn, rn.child[ANN_HI]);
  END_NEW_TIME(19)

  // Sort closer branches
  ADD_BENCHMARK_DESCRIPTION("Storing and sorting node scores", 17)
  START_TIME()
  NodeScore scores[2] = { NodeScore(&qn, rn.child[ANN_LO], s_l), NodeScore(&qn, rn.child[ANN_HI], s_r) } ;
  sort2_sn_stable(scores); // Sort by score using sorting network
  END_TIME(17)

  // Prioritize closer branches, return at first infinity. Always update the bounds after the DFS.
  const int n_branches = 2;
  for (int i = 0; i < n_branches; ++i){
    if (scores[i].score != ANN_DIST_INF){
      DFS(*scores[i].lhs, *scores[i].rhs);
      ADD_BENCHMARK_DESCRIPTION("Updating split node bounds", 13)
      START_TIME()
      R_OUT("Updating bounds--> q:%d, r:%d\n", qn.id, rn.id);
      UpdateBounds(&qn, &rn);
      END_TIME(13)
    } else {
      break; // break out of loop and end recursion at first infinity
    }
  }
  return;
  // int cond = ((scores[0].score == ANN_DIST_INF ? COND_0 : 0) | (scores[1].score == ANN_DIST_INF ? COND_1 : 0));
  // switch(cond){
  //   case 0:
  //     PathDFS(lpath, rpath, scores[0]), UpdateBounds(&rn);
  //     // R_OUT("-- Continuing from parent: (q:%d, r:%d) --> q_: %d, r_: %d\n", getNodeID(&qn), getNodeID(&rn), getNodeID(scores[1].lhs), getNodeID(scores[1].rhs));
  //     PathDFS(lpath, rpath, scores[1]), UpdateBounds(&rn);
  //     break;
  //   case COND_1:
  //     R_OUT("Cond 1: Recurse 0 (%d, %d), Prune 1 (%d, %d)\n", scores[0].lhs, scores[0].rhs, scores[1].lhs, scores[1].rhs);
  //     PathDFS(lpath, rpath, scores[0]), UpdateBounds(&rn);
  //     break;
  //   case COND_0 + COND_1: return; break;
  // }
  // UpdateBounds(&rn); // Only need to update split node
};

template <class METRIC_T>
void DualTreeSearch<METRIC_T>::DFS(ANNkd_leaf& qn, ANNkd_leaf& rn){
  R_OUT("leaf DFS: (q:%d, r:%d)\n", getQN_ID(&qn), getRN_ID(&rn));
  // auto key = &qn < &rn ? std::make_pair(&qn, &rn) : std::make_pair(&rn, &qn);
  ADD_BENCHMARK_DESCRIPTION("Searching cache for key", 14)
  START_TIME()
  auto key = std::make_pair(&qn, &rn);
  CACHE::iterator cache_hit = cache_map.find(key);
  END_TIME(14)
  if (cache_hit == cache_map.end()){ // node combination has never been visited
    R_OUT("-- VISITING Regular Base case -- %d\n", 0);
    BaseCase(&qn, &rn); // Base case updates bounds and cache for leaf nodes
  }
  else {
    // This pair has been visited before, determine whether the distance map was saved
    // in the current (qn, rn) order, or if it was saved under the (rn, qn) order.
    R_OUT("-- VISITING CACHE (ordered: %d)-- \n", (*cache_hit).second.reverse_ordered);
    BaseCaseCached(&qn, &rn);
  }
}

template <class METRIC_T>
ANNdist DualTreeSearch<METRIC_T>::computeBoxDist(ANNkd_node* qn, ANNkd_node* rn){
  ANNdist dist = 0., d, t1, t2;
  for (int d_i = 0; d_i < qtree->dim; ++d_i){
    t1 = rinfo[rn->id].lo[d_i] - qinfo[qn->id].hi[d_i];
    t2 = qinfo[qn->id].lo[d_i] - rinfo[rn->id].hi[d_i];
    d = ((t1 + std::abs(t1)) + (t2 + std::abs(t2)))*0.5;
    dist += ANN_POW(d);
    // R_OUT("t1: %f, t2: %f, box dist: %f\n", t1, t2, dist);
    // switch(cut_dir[d_i]){
    //   case CUT_LO:
    //     t = node_info[rn].lo[d_i] - node_info[qn].hi[d_i];
    //     dist += (t * t);
    //     // metric.distance(node_info[rn].lo, node_info[qn].hi, d_i, dist);
    //     break;
    //   case CUT_HI:
    //     t = node_info[qn].lo[d_i] - node_info[rn].hi[d_i];
    //     dist += (t * t);
    //     // metric.distance(node_info[qn].lo, node_info[rn].hi, d_i, dist);
    //     break;
    //   case CUT_OVERLAP: // don't include contributions of overlapping boxes
    //     break;
    // }
  }
  R_OUT("Box distance: (%d, %d) --> %f, cut dir: [%d, %d]\n", getQN_ID(qn), getRN_ID(rn), dist, cut_dir[0], cut_dir[1]);
  // return(comp > 1 ? sqrt(dist) : dist);
  return(dist);
}


template <class METRIC_T>
List DualTreeSearch<METRIC_T>::getTreeProperties(bool query){
  // Retrieve all properties built from the tree as a List
  List tree_properties = List();
  NODE_INFO ni_info = query ? qinfo : rinfo;
  for(NODE_INFO::iterator ni = ni_info.begin(); ni != ni_info.end(); ++ni) {
    // bool is_leaf = dynamic_cast<ANNkd_leaf*>(kv.first) != NULL;
    List tree_prop = List::create(_["id"] = (*ni).id,
                                  _["bound"] = (*ni).bound,
                                  _["centroid"] = ptToVec((*ni).centroid, qtree->dim),
                                  _["radius"] = (*ni).max_radius,
                                  _["lb"] = ptToVec((*ni).lo, qtree->dim), _["ub"] =  ptToVec((*ni).hi, qtree->dim));
    // if (is_leaf){
    //   ANNkd_leaf* leaf_node = dynamic_cast<ANNkd_leaf*>(kv.first);
    //   tree_prop["idx"] = idxArrayToVec(leaf_node->bkt, leaf_node->n_pts);
    // }
    tree_properties.push_back(tree_prop);
  }
  return(tree_properties);
}
// if (q_hi == 19 && r_lo == 31){
//   R_OUT("****scores: %f, %f, %f, %f\n", scores[0].score, scores[1].score, scores[2].score, scores[3].score);
// }
//
// // Prioritize closer branches, return at first infinity. Always update the bounds after the DFS.
// int cond = ((scores[0].score == ANN_DIST_INF ? COND_0 : 0) | (scores[1].score == ANN_DIST_INF ? COND_1 : 0) |
//             (scores[2].score == ANN_DIST_INF ? COND_2 : 0) | (scores[3].score == ANN_DIST_INF ? COND_3 : 0));
// if (q_hi == 19 && r_lo == 31){
//   R_OUT("****Condition: %d\n", cond);
// }
// switch(cond){
// case 0:
//   R_OUT("Base cond: Recurse all\n", 0);
//   PathDFS(lpath, rpath, scores[0]), UpdateBounds(&qn, &rn);
//   // R_OUT("-- Continuing from parent: (q:%d, r:%d) --> q_: %d, r_: %d\n", getNodeID(&qn), getNodeID(&rn), getNodeID(scores[1].lhs), getNodeID(scores[1].rhs));
//   PathDFS(lpath, rpath, scores[1]), UpdateBounds(&qn, &rn);
//   // R_OUT("-- Continuing from parent: (q:%d, r:%d) --> q_: %d, r_: %d\n", getNodeID(&qn), getNodeID(&rn), getNodeID(scores[2].lhs), getNodeID(scores[2].rhs));
//   PathDFS(lpath, rpath, scores[2]), UpdateBounds(&qn, &rn);
//   // R_OUT("-- Continuing from parent: (q:%d, r:%d) --> q_: %d, r_: %d\n", getNodeID(&qn), getNodeID(&rn), getNodeID(scores[3].lhs), getNodeID(scores[3].rhs));
//   PathDFS(lpath, rpath, scores[3]), UpdateBounds(&qn, &rn);
//   break;
// case COND_3:
//   R_OUT("Condition 3: Recurse 0, 1, 2\n", 0);
//   PathDFS(lpath, rpath, scores[0]), UpdateBounds(&qn, &rn);
//   // R_OUT("-- Continuing from parent: (q:%d, r:%d) --> q_: %d, r_: %d\n", getNodeID(&qn), getNodeID(&rn), getNodeID(scores[1].lhs), getNodeID(scores[1].rhs));
//   PathDFS(lpath, rpath, scores[1]), UpdateBounds(&qn, &rn);
//   // R_OUT("-- Continuing from parent: (q:%d, r:%d) --> q_: %d, r_: %d\n", getNodeID(&qn), getNodeID(&rn), getNodeID(scores[2].lhs), getNodeID(scores[2].rhs));
//   PathDFS(lpath, rpath, scores[2]), UpdateBounds(&qn, &rn);
//   break;
// case COND_2 + COND_3:
//   R_OUT("Condition 2+3: Recurse 0, 1\n", 0);
//   PathDFS(lpath, rpath, scores[0]), UpdateBounds(&qn, &rn);
//   // R_OUT("-- Continuing from parent: (q:%d, r:%d) --> q_: %d, r_: %d\n", getNodeID(&qn), getNodeID(&rn), getNodeID(scores[1].lhs), getNodeID(scores[1].rhs));
//   PathDFS(lpath, rpath, scores[1]), UpdateBounds(&qn, &rn);
//   break;
// case COND_1 + COND_2 + COND_3:
//   R_OUT("Condition 1+2+3: Recurse 0\n", 0);
//   PathDFS(lpath, rpath, scores[0]), UpdateBounds(&qn, &rn);
// case COND_0 + COND_1 + COND_2 + COND_3:
//   R_OUT("Condition PRUNE ALL\n", 0);
//   return;
//   break;
// }
// UpdateBounds(&qn, &rn); // update both nodes
// }

// DualTreeSearch<L_2> a; // (nullptr, nullptr, NODE_INFO(), L_2(0), 0);
template class DualTreeSearch<L_2>;
template class DualTreeSearch<L_1>;
template class DualTreeSearch<L_inf>;
template class DualTreeSearch<L_p>;
template class DualTreeSearch<RSL>;
