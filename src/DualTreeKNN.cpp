#include "DualTreeKNN.h"

// template <class METRIC_T>
// DualTreeKNN<METRIC_T>::DualTreeKNN() : DualTreeSearch<METRIC_T>::DualTreeSearch(), k(0) {};

template <class METRIC_T>
DualTreeKNN<METRIC_T>::DualTreeKNN(ANNkd_tree* qtree, ANNkd_tree* rtree, const NODE_INFO& qinfo, const NODE_INFO& rinfo, METRIC_T& m, const int k)
  : DualTreeSearch<METRIC_T>::DualTreeSearch(qtree, rtree, qinfo, rinfo, m), k(k)  {
  knn_pq = std::vector<ANNmin_k*>(qtree->n_pts); // priority queue for the closest k points
  for (int i = 0; i < qtree->n_pts; ++i){ knn_pq[i] = new ANNmin_k(k); } // initialize pr queue

  // Benchmarking / Profiling
  benchmarks = std::vector<double>(20, 0.0);
  benchmark_descrip = std::vector<std::string>(20);
  dts::n_comparisons = 0;
  dts::n_pruned = 0;
}

// Update function for split nodes
template <class METRIC_T>
void DualTreeKNN<METRIC_T>::UpdateBounds(ANNkd_split* node, NODE_INFO& ninfo){
  START_TIME()
  double left_bound = ninfo[node->child[ANN_LO]->id].bound;
  double right_bound = ninfo[node->child[ANN_HI]->id].bound;
  ninfo[node->id].bound = std::max(left_bound, right_bound);
  END_TIME(7)
  return;
}

template <class METRIC_T>
void DualTreeKNN<METRIC_T>::UpdateBounds(ANNkd_leaf* node, NODE_INFO& ninfo) { return; } // implemented in base case


  // double max_dist = 0;
  // for (int i = 0; i < node->n_pts; ++i){
  //   int idx = node->bkt[i];
  //   double max_knn_dist = knn_pq[idx]->max_key(); // max knn dist
  //   if (max_knn_dist > max_dist){ max_dist = max_knn_dist; }
  // }
  // if (max_dist > node_info[node].bound) { node_info[node].bound = max_dist; }
  // Rcout << "Updating leaf bound\n";
  // for (int i = 0; i < node->n_pts; ++i){
  //   knn_pq[node->bkt[i]]->max_key();
  // }
  // std::unordered_map<ANNkd_node*, std::vector<ANNmin_k::mk_node*>>::const_iterator knn_dist = knn_ptr.find(node);
  // if (knn_dist != knn_ptr.end()){
  //   ANNmin_k::mk_node* max_node = *std::max_element(knn_dist->second.begin(), knn_dist->second.end(), node_key);
  //   node_info[node].bound = max_node->key;
  //   for (int i = 0; i < node->n_pts; ++i){
  //     R_OUT("true ptr: %p\n", (void*) knn_pq[node->bkt[i]]);
  //     R_OUT("true ptr: %p\n", (void*) &knn_pq[node->bkt[i]]->mk[k-1]);
  //     R_OUT("true knn_dist: %f\n", knn_pq[node->bkt[i]]->max_key());
  //   }
  //   for (int i = 0; i < knn_ptr.at(node).size(); ++i){
  //     R_OUT("my ptr: %p\n", (void*) knn_ptr.at(node)[i]);
  //     R_OUT("knn_dist: %f\n", knn_ptr.at(node)[i]->key);
  //   }
  //   R_OUT("Node updated bound: %f\n", max_node->key);
  // }
// }

// KNN Base Case: Compute all the pairwise distances for each pair of points between the leaf the nodes
template <class METRIC_T>
void DualTreeKNN<METRIC_T>::BaseCase(ANNkd_leaf* qn, ANNkd_leaf* rn){
  const bool same_node = qn == rn;
  double d_pq;
  double q_max_el = 0; // The maximum KNN distances for the query and reference nodes

  // Number of pointwise comparisons about to make
  const int n_pt_comparisons = same_node ? (((qn->n_pts * qn->n_pts + 1) / 2) - qn->n_pts) : qn->n_pts * rn->n_pts;
  dts::n_comparisons += n_pt_comparisons;

  // Create a cache for the distances
  std::vector<ANNdist>* dist_cache = new std::vector<ANNdist>();
  dist_cache->reserve(n_pt_comparisons);

  START_TIME()
  for (int i = 0; i < qn->n_pts; ++i){
    R_OUT("Starting base case comparisons for q pt index %d (n_pts = %d)\n", i, qn->n_pts);
    int q_idx = qn->bkt[i];
    if (q_idx > knn_pq.size()){
      R_OUT("Trying to access index %d from a vector of size %d\n", q_idx, (int) knn_pq.size());
      stop("here3");
    }
    d_pq = knn_pq[q_idx]->max_key();
    for (int j = (same_node ? i+1 : 0); j < rn->n_pts; ++j){ // If same node, all n x n pairwise distances would be redundant
      { R_OUT("Starting base case comparisons for r pt index %d (n_pts = %d)\n", j, rn->n_pts); }
      int r_idx = rn->bkt[j];

      // TODO: change when qtree != rtree
      //if (q_idx == r_idx) continue;
      // double dist = dts::computeDist(dts::qtree->pts[q_idx], dts::rtree->pts[r_idx], dts::qtree->dim);
      if (q_idx > dts::qtree->n_pts || r_idx > dts::rtree->n_pts){ stop("Q or R idx greater than number of points in tree"); }
      ANNdist dist = dts::metric.distance(q_idx, r_idx);
      R_OUT("Comparing pts: (q:%d, r:%d, dist:%f, q_max_knn: %f)\n", q_idx, r_idx, dist, knn_pq.at(q_idx)->max_key());

      if (dist < d_pq){ knn_pq.at(q_idx)->insert(dist, r_idx); } // if the distance is low enough, insert into priority queue for query pt
      if (same_node && dist < knn_pq.at(r_idx)->max_key()){ knn_pq.at(r_idx)->insert(dist, q_idx); } // exception case when at the same node
      dist_cache->push_back(dist); // Store the newly computed things in the cache
    }
    if (knn_pq.at(q_idx)->max_key() > q_max_el){
      q_max_el = knn_pq.at(q_idx)->max_key();
    }
  }
  END_TIME(0)

  // Update the query node bound w/ the max KNN dist
  dts::setQN_Bound(qn, q_max_el);

  START_TIME()
  dts::cache_map.emplace(std::make_pair(qn, rn), cache(false, dist_cache)); // emplace a new cache object
  dts::cache_map.emplace(std::make_pair(rn, qn), cache(true, dist_cache)); // emplace a new cache object, order reversed
  END_TIME(3)
}

template <class METRIC_T>
void DualTreeKNN<METRIC_T>::BaseCaseCached(ANNkd_leaf* qn, ANNkd_leaf* rn){
  double d_pq;
  double q_max_el = 0; // The maximum KNN distances for the query and reference nodes
  START_TIME()
  auto key = std::make_pair(qn, rn);
  const std::vector<ANNdist>& leaf_dist = *dts::cache_map.at(key).dist_map;
  const bool rev_ordered = dts::cache_map.at(key).reverse_ordered;
  END_TIME(2)

  dts::n_accesses += qn->n_pts * rn->n_pts;

  START_TIME()
  int counter = 0;
  for (int i = 0; i < qn->n_pts; ++i){
    int q_idx = qn->bkt[i]; // get query point
    for (int j = 0; j < rn->n_pts; ++j){
      int r_idx = rn->bkt[j]; // get the reference point
      ANNdist dist = leaf_dist.at(rev_ordered ? j*qn->n_pts + i : counter++);
      R_OUT("Dist between pts: (q:%d, r:%d) --> %f\n", q_idx, r_idx, dist);
      d_pq = knn_pq[q_idx]->max_key(); // get the current knn distance
      if (dist < d_pq){ knn_pq[q_idx]->insert(dist, r_idx); }  // if the distance is low enough, insert into priority queue for query pt
    }

    // Update bound
    if (knn_pq[q_idx]->max_key() > dts::getQN_Bound(qn)){
      dts::setQN_Bound(qn, knn_pq[q_idx]->max_key()); // update max knn; don't use d_pq as might be out of date
    }
  }
  END_TIME(1)
}


// KNN Score Function
template <class METRIC_T>
ANNdist DualTreeKNN<METRIC_T>::Score(ANNkd_node* qn, ANNkd_node* rn){
  if (qn == rn){ return(0.0); } // trivial case: nodes that are the same have 0 distance to each other
  double q_bnd = dts::getQN_Bound(qn); // If doing regular dual traversal (std::max(getBound(qn), getBound(rn)); // experimental triple traversal)

  // Compute the bound based on the metric
  ANNdist d_min = dts::computeBound(qn, rn);
  R_OUT("Comparing: [%d, %d] (box_lb: %f, current bound: %f) \n", dts::getQN_ID(qn), dts::getRN_ID(rn), d_min, q_bnd);
  if(d_min > q_bnd) {
    dts::n_pruned++;
    R_OUT("... pruned [%d, %d] (box_lb: %f, current bound: %f) \n", dts::getQN_ID(qn), dts::getRN_ID(rn), d_min, q_bnd);
    return(ANN_DIST_INF); // distance between bounding boxes of nodes is high enough: prune combination
  }
  return(d_min);

  // Using the circumsphere bound
  // ANNpoint q_centroid = node_info[qn].centroid; // getCentroid(qn);
  // ANNpoint r_centroid = node_info[rn].centroid; //getCentroid(rn)
  // R_OUT("Centroid dist: %f, max qn: %f, max rn: %f\n", annDist(qtree->dim, q_centroid, r_centroid), std::pow(getMaxRadius(qn), 2), std::pow(getMaxRadius(rn), 2));
  // double d_min = std::sqrt(annDist(qtree->dim, q_centroid, r_centroid)) - getMaxRadius(qn) - getMaxRadius(rn);
  // d_min = d_min < 0 ? 0 : std::pow(d_min, 2);
  // if (d_min > q_bnd){
  //   n_pruned++;
  //   // Rprintf("... pruned [%d, %d] (box_lb: %f, current bound: %f) \n", getNodeID(qn), getNodeID(rn), d_min, getBound(qn));
  //   return(ANN_DIST_INF); // distance between bounding boxes of nodes is high enough: prune combination
  // }
  // return(d_min);
  // R_OUT("(%d, %d) Circumsphere dist: %f, True Box dist: %f, Q Bound: %f\n", getNodeID(qn), getNodeID(rn), d_min, actual_box_dist, q_bnd);

}


// Convert the results in the priority queue suitable for returning over to R (1-based)
template <class METRIC_T>
List DualTreeKNN<METRIC_T>::getKnnResults(){
  const int n = dts::qtree->n_pts;
  IntegerMatrix idx = Rcpp::no_init_matrix(n, k);
  NumericMatrix dist = Rcpp::no_init_matrix(n, k);
  for (int i = 0; i < n; ++i){
    for (int k_i = 0; k_i < k; ++k_i){
      dist(i, k_i) = (ANNdist) knn_pq[i]->ith_smallest_key(k_i);
      idx(i, k_i) = (int) knn_pq[i]->ith_smallest_info(k_i);
    }
    // if (i == 0){
    //   Rcout << "0 pt nearest neighbor: " << dist(i, 0) << std::endl;
    // }
  }
  return(List::create(_["dist"] = dist, _["idx"] = idx + 1));
}

// Explicit instantiations
template class DualTreeKNN<L_2>;
template class DualTreeKNN<L_1>;
template class DualTreeKNN<L_inf>;
template class DualTreeKNN<L_p>;
template class DualTreeKNN<RSL>;
// DualTreeKNN<L_2> a;
