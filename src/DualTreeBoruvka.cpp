#include "DualTreeBoruvka.h"

std::vector<double> benchmarks;
std::vector<std::string> benchmark_descrip;
std::chrono::time_point<std::chrono::high_resolution_clock> start;
std::chrono::time_point<std::chrono::high_resolution_clock> finish;

template <class METRIC_T>
DualTreeBoruvka<METRIC_T>::DualTreeBoruvka(ANNkd_tree* qtree, ANNkd_tree* rtree, const NODE_INFO& qinfo, const NODE_INFO& rinfo, METRIC_T& m)
  : DualTreeSearch<METRIC_T>::DualTreeSearch(qtree, rtree, qinfo, rinfo, m), CC(UnionFind(qtree->n_pts)) {
  N = std::vector<double_edge>();
  N.reserve(qtree->n_pts);
  for (int i = 0; i < qtree->n_pts; ++i){ N.push_back(double_edge(i, i, ANN_DIST_INF)); }
  ALL_CC_SAME = std::vector<int>(dts::qinfo.size(), -1); // By default, every node (and descendents) are in different components
  // for(auto& kv : dts::node_info) { ALL_CC_SAME[kv.first] = -1; } // By default, every node (and descendents) are in different components
  benchmarks = std::vector<double>(20, 0.0);
  benchmark_descrip = std::vector<std::string>(20);
  dts::n_comparisons = 0;
  dts::n_pruned = 0;
  n_accesses = 0;
  // dist_map = std::unordered_map< pair<int,int>, ANNdist >();
}

template <class METRIC_T>
void DualTreeBoruvka<METRIC_T>::BaseCase(ANNkd_leaf* qn, ANNkd_leaf* rn){
  const bool same_node = qn == rn;

  // Number of pointwise comparisons about to make
  const int n_pt_comparisons = same_node ? (((qn->n_pts * qn->n_pts + 1) / 2) - qn->n_pts) : qn->n_pts * rn->n_pts;
  dts::n_comparisons += n_pt_comparisons;

  // Create a cache for the distances
  ADD_BENCHMARK_DESCRIPTION("Reserving distance cache", 0)
  START_TIME()
  std::vector<ANNdist>* dist_cache = new std::vector<ANNdist>();
  dist_cache->reserve(n_pt_comparisons);
  END_TIME(0)

  // To allow the DTB to score nodes, it's necessary to keep track of whether all points fall within the same component for a given node
  bool q_same_comp = true, r_same_comp = true; // check to see if all points in both nodes belong to same component
  const int q_orig_comp = CC.Find(qn->bkt[0]), r_orig_comp = CC.Find(rn->bkt[0]);

  // To update the query nodes bound appropriately, need effectively the max(D[pq]_knn) forall q in Q, where k = 1.
  // Rather than instantiate a priroity queue, simply keep track of the smallest pairwise distance for each query
  // point, and then take the maximum of these to get the bound.
  double qpt_nearest_dist = ANN_DIST_INF, qnode_max_dist = -1;
  for (int i = 0; i < qn->n_pts; ++i){
    int q_idx = qn->bkt[i]; // query point index
    const int q_comp = CC.Find(q_idx); // query point's component
    for (int j = (same_node ? i+1 : 0); j < rn->n_pts; ++j){
      int r_idx = rn->bkt[j]; // reference point index
      const int r_comp = CC.Find(r_idx); // Reference points component

      // Compute distance
      ADD_BENCHMARK_DESCRIPTION("Computing metric distance", 4)
      START_TIME()
      ANNdist dist = dts::metric.distance(q_idx, r_idx);
      R_OUT("Comparing pts: (q_pt:%d, r_pt:%d) --> %f\n", q_idx, r_idx, dist);
      dist_cache->push_back(dist); // Store the newly computed things in the cache
      END_TIME(4)


      // TODO: change when qtree != rtree
      // If query and refernce points lie in different connected components and their distance is
      // less than the current query components lowest distance (N[q_comp].weight), then update the edge.
      if (q_comp != r_comp && dist < N[q_comp].weight){
        // R_OUT("Pt query comp: %d, pt ref comp: %d\n", CC.Find(q_idx), CC.Find(r_idx));
        R_OUT("Update Q Component F([%d, %d]) shortest edge: (pt: %d, pt: %d, %f)\n", q_comp, r_comp, q_idx, r_idx, dist);
        ADD_BENCHMARK_DESCRIPTION("Updating Q components shortest edge", 5)
        START_TIME();
        N[q_comp] = { q_idx, r_idx, dist };
        END_TIME(5);
      }

      // If the query and reference nodes are the same, and the query and reference components
      // are different, update the reference components nearest edge as well.
      // TODO: change when qtree != rtree
      if (same_node){
        ADD_BENCHMARK_DESCRIPTION("Updating R components shortest edge", 15)
        START_TIME();
        double d_pr = N[r_comp].weight;  // get current lowest weight edge
        if (dist < d_pr && q_comp != r_comp){
          R_OUT("Update R Component F([%d, %d]) shortest edge: (%d, %d, %f)\n", r_comp, q_comp, q_idx, r_idx, dist);
          N[r_comp] = { r_idx, q_idx, dist };
          d_pr = dist;
        }
        END_TIME(15);
      }

      // For the first iteration only, check to see if every reference point lies in the same connected component.
      if (i == 0 && r_same_comp){ r_same_comp = (r_orig_comp == r_comp); }

      // Update nearest neighbor distance for given point
      // IMPORTANT: Since this distance ends up possibly being the bound of the query node,
      // must only lower this distance if query and reference points are in separate components!
      if (q_comp != r_comp && dist < qpt_nearest_dist){ qpt_nearest_dist = dist; }
    }
    ADD_BENCHMARK_DESCRIPTION("Checking component conditions", 7)
    START_TIME();
    // Check to see if every query point lies in the same component
    if (q_same_comp){ q_same_comp = (q_orig_comp == q_comp); }

    // Update the max nearest neighbor distance for the given point
    if (qpt_nearest_dist > qnode_max_dist){ qnode_max_dist = qpt_nearest_dist; }
    qpt_nearest_dist = ANN_DIST_INF; // reset
    END_TIME(7);
  }


  ADD_BENCHMARK_DESCRIPTION("Creating caches, updating leaf node bounds, updating component status", 8)
  START_TIME();
  // Update the final status on whether all points are part of same component
  ALL_CC_SAME[qn->id] = q_same_comp ? q_orig_comp : -1;
  ALL_CC_SAME[rn->id] = r_same_comp ? r_orig_comp : -1;

  // Update query node bound with the maximum k=1 nearest neighbor distance of all query points
  dts::qinfo[qn->id].bound = (qnode_max_dist == -1) ? ANN_DIST_INF : qnode_max_dist; // If no neighbors were found, use INF

  // Create the caches
  dts::cache_map.emplace(std::make_pair(qn, rn), cache(false, dist_cache)); // emplace a new cache object, regular ordering
  dts::cache_map.emplace(std::make_pair(rn, qn), cache(true, dist_cache)); // emplace a new cache object, order reversed
  END_TIME(8);
}

template <class METRIC_T>
void DualTreeBoruvka<METRIC_T>::BaseCaseCached(ANNkd_leaf* qn, ANNkd_leaf* rn){
  // R_OUT("Checking components: [%d, %d], comp: %d, r0 comp: %d\n", q_orig_comp, r_orig_comp);

  ADD_BENCHMARK_DESCRIPTION("Making key pairs", 2)
  START_TIME()
  auto key = std::make_pair(qn, rn);
  const std::vector<ANNdist>& leaf_dist = *dts::cache_map.at(key).dist_map;
  const bool rev_ordered = dts::cache_map.at(key).reverse_ordered;
  END_TIME(2)

  const bool same_node = qn == rn;
  int leaf_c = 0;
  n_accesses += qn->n_pts * rn->n_pts; // update number of cache accesses

  // To allow the DTB to score nodes, it's necessary to keep track of whether all points fall within the same component for a given node
  bool q_same_comp = true, r_same_comp = true; // check to see if all points in both nodes belong to same component
  const int q_orig_comp = CC.Find(qn->bkt[0]), r_orig_comp = CC.Find(rn->bkt[0]);

  // To update the query nodes bound appropriately, need effectively the max(D[pq]_knn) forall q in Q, where k = 1.
  // Rather than instantiate a priroity queue, simply keep track of the smallest pairwise distance for each query
  // point, and then take the maximum of these to get the bound.
  double qpt_nearest_dist = ANN_DIST_INF, qnode_max_dist = -1;

  ADD_BENCHMARK_DESCRIPTION("Doing base case, cached", 1)
  START_TIME()
  int counter = 0;
  for (int i = 0; i < (same_node ? qn->n_pts - 1 : qn->n_pts); ++i){
    int q_idx = qn->bkt[i]; // get query point
    const int q_comp = CC.Find(q_idx); // get the query component
    for (int j = (same_node ? i+1 : 0); j < rn->n_pts; ++j){
      int r_idx = rn->bkt[j]; // get the reference point
      const int r_comp = CC.Find(r_idx); // get the reference component
      ANNdist dist = leaf_dist.at(rev_ordered ? j*qn->n_pts + i : counter++); // get cached distance; note the ordering
      R_OUT("Comparing pts: (q_pt:%d, r_pt:%d) --> %f\n", q_idx, r_idx, dist);
      if (r_comp != q_comp && dist < N[q_comp].weight){
        R_OUT("Update Q Component F([%d, %d]) shortest edge: (%d, %d, %f)\n", q_comp, r_comp, q_idx, r_idx, dist);
        N[q_comp] = { q_idx, r_idx, dist };
      }
      if (same_node && (q_comp != r_comp && dist < N[r_comp].weight)){
        R_OUT("Update R Component F([%d, %d]) shortest edge: (%d, %d, %f)\n", r_comp, q_comp, q_idx, r_idx, dist);
        N[r_comp] = { r_idx, q_idx, dist };
      }

      // For the first iteration only, check to see if every reference point lies in the same connected component.
      if (i == 0 && r_same_comp){ // R_OUT("R same comp: %d == %d\n", r_orig_comp, r_comp);
        r_same_comp = (r_orig_comp == r_comp);
      }
      // Update nearest neighbor distance for given point
      if (q_comp != r_comp && dist < qpt_nearest_dist){
        R_OUT("Updating q nearest dist qp_nearest: %f --> %f\n", qpt_nearest_dist, dist);
        qpt_nearest_dist = dist;
      }
    }

    // Check to see if every query point lies in the same component
    if (q_same_comp){
      // R_OUT("Q same comp: %d == %d\n", q_orig_comp, q_comp);
      q_same_comp = (q_orig_comp == q_comp);
    }

    // Update the max nearest neighbor distance for the given point
    // R_OUT("Qpt:%d nearest dist: %f, qnode max dist: %f\n", q_idx, qpt_nearest_dist, qnode_max_dist);
    if (qpt_nearest_dist > qnode_max_dist){ qnode_max_dist = qpt_nearest_dist; }
    qpt_nearest_dist = ANN_DIST_INF; // reset
  }
  END_TIME(1)

  // Update the final status on whether all points are part of same component
  R_OUT("Component status: Q[%d]:%d->%d, R[%d]:%d->%d, (q max k=1 nn dist: %f)\n",
        dts::getQN_ID(qn), ALL_CC_SAME[qn->id], q_same_comp ? q_orig_comp : -1,
        dts::getRN_ID(rn), ALL_CC_SAME[rn->id], r_same_comp ? r_orig_comp : -1,
        qnode_max_dist);

  ADD_BENCHMARK_DESCRIPTION("Updating component status, cached", 16)
  START_TIME()
  ALL_CC_SAME[qn->id] = q_same_comp ? q_orig_comp : -1;
  ALL_CC_SAME[rn->id] = r_same_comp ? r_orig_comp : -1;

  // Update query node bound with the maximum k=1 nearest neighbor distance of all query points
  dts::qinfo[qn->id].bound = (qnode_max_dist == -1) ? ANN_DIST_INF : qnode_max_dist; // If no neighbors were found, use INF
  END_TIME(16)
  return;
}


template <class METRIC_T>
ANNdist DualTreeBoruvka<METRIC_T>::Score(ANNkd_node* qn, ANNkd_node* rn){
  if (qn == rn){ return(0); } // trivial case: nodes that are the same have 0 distance to each other
  double q_bnd = dts::getQN_Bound(qn);

  // Compute minimum distance between qn and rn
  ANNdist d_min = dts::computeBound(qn, rn);

  if (d_min < q_bnd){
    const int q_comp = ALL_CC_SAME[qn->id], r_comp = ALL_CC_SAME[rn->id];
    // R_OUT("Min dist between nodes %d and %d: %f, Q comp: %d, R comp: %d\n", node_info[qn].id, node_info[rn].id, d_min, q_comp, r_comp);
    if (q_comp == r_comp && q_comp != -1){
      dts::n_pruned++;
      R_OUT("* Pruning: (id: %d, comp: %d), (id: %d, comp: %d)\n", dts::qinfo[qn->id].id, q_comp, dts::rinfo[rn->id].id, r_comp);
      return(ANN_DIST_INF); // prune since all descendent points in both nodes are in the same component
    }
    return(d_min); // o/w return distance between convex subsets to prioritize recursion
  }
  R_OUT("* Pruning: (id: %d, comp: %d), (id: %d, comp: %d), { box_dist: %f, q_bnd: %f, res: %d}\n",
        dts::qinfo[qn->id].id, ALL_CC_SAME[qn->id], dts::rinfo[rn->id].id, ALL_CC_SAME[rn->id],
                                                                                  d_min, q_bnd, d_min < q_bnd);
  dts::n_pruned++;
  return(ANN_DIST_INF); // prune since minimum distance between nodes is greater than bound
} // TODO: Update Bounds after every DFS, and add ALL_CC_SAME to bound update!

template <class METRIC_T>
void DualTreeBoruvka<METRIC_T>::UpdateBounds(ANNkd_split* node, NODE_INFO& ninfo){
  // double left_bound = dts::getBound(node->child[ANN_LO]);
  // double right_bound = dts::getBound(node->child[ANN_HI]);
  // double new_bound = std::max(left_bound, right_bound);
  // if (new_bound > dts::node_info[node].bound){
  //   dts::node_info[node].bound = new_bound;
  // }
  ninfo[node->id].bound = std::max(ninfo[node->child[ANN_LO]->id].bound, ninfo[node->child[ANN_HI]->id].bound);
  if (ALL_CC_SAME[node->child[ANN_LO]->id] != -1){
    ALL_CC_SAME[node->id] = (ALL_CC_SAME[node->child[ANN_LO]->id] == ALL_CC_SAME[node->child[ANN_HI]->id]) ? ALL_CC_SAME[node->child[ANN_LO]->id] : -1;
  }
}

// Unimplemented; leaf bound updates handled in base case
template <class METRIC_T>
void DualTreeBoruvka<METRIC_T>::UpdateBounds(ANNkd_leaf* node, NODE_INFO& ninfo){}


// To actually compute the minimum spanning tree, need to run DFS multiple times
template <class METRIC_T>
NumericMatrix DualTreeBoruvka<METRIC_T>::MST(){
  // Calculate the minimum spanning tree
  const int n_edges = dts::qtree->n_pts - 1;
  NumericMatrix mst = Rcpp::no_init_matrix(n_edges, 3);
  bool is_fully_connected = false;
  int i = 0, n_iter = 1;
  do {
    Rcpp::checkUserInterrupt();

    ADD_BENCHMARK_DESCRIPTION("Resetting node bounds", 9)
    START_TIME()
    dts::resetBounds(); // Must reset the bounds between each boruvka step
    END_TIME(9)

    // Depth-First Traversal solves the Boruvka step
    ADD_BENCHMARK_DESCRIPTION("Doing Boruvka Steps", 18)
    START_NEW_TIME()
    dts::DFS();
    END_NEW_TIME(18)

    // Union the components together
    ADD_BENCHMARK_DESCRIPTION("Sorting edges", 10)
    START_TIME()
    std::sort(N.begin(), N.end()); // sort nearest-component edges in increasing order
    END_TIME(10)

    R_OUT("Looking at edges\n", "");
    ADD_BENCHMARK_DESCRIPTION("Unioning edges", 3)
    START_TIME()
    for (std::vector<double_edge>::iterator e = N.begin(); e != N.end(); ++e){
      const double_edge& c_edge = (*e);
      if (CC.Find(c_edge.from) != CC.Find(c_edge.to) && c_edge.weight != ANN_DIST_INF){
        R_OUT("Unioning: (%d, %d, %f)\n", c_edge.from, c_edge.to, c_edge.weight);
        CC.Union(c_edge.from, c_edge.to); // union the components
        mst(i++, _) = NumericVector::create(c_edge.from, c_edge.to, c_edge.weight); // add the edge to the spanning tree
      }
      (*e).weight = ANN_DIST_INF; // Reset edge weight
    }
    END_TIME(3)

    ADD_BENCHMARK_DESCRIPTION("Checking if fully connected", 6)
    START_TIME()
    is_fully_connected = fully_connected();
    END_TIME(6)
    n_iter++;
  } while (!is_fully_connected); // continue until we have a spanning tree

  // Return the final spanning tree
  return(mst);
}

template <class METRIC_T>
NumericMatrix DualTreeBoruvka<METRIC_T>::getResults(){
  const int n = dts::qtree->n_pts;
  NumericMatrix mst = Rcpp::no_init_matrix(n, 3);
  int i = 0;
  for (std::vector<double_edge>::iterator edge_i = N.begin(); edge_i != N.end(); ++edge_i){
    // if ((*edge_i).weight != ANN_DIST_INF){
      mst(i++, _) = NumericVector::create((*edge_i).from, (*edge_i).to, (*edge_i).weight);
    // }
  }
  return(mst);
}

// explicit instantiations
template class DualTreeBoruvka<L_2>;
template class DualTreeBoruvka<L_1>;
template class DualTreeBoruvka<L_inf>;
template class DualTreeBoruvka<L_p>;
template class DualTreeBoruvka<RSL>;
