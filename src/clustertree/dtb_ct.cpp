#include "dtb_ct.h" // Dual Tree Borvuka definitions

// use default constructor
DTB_CT::DTB_CT(const bool prune, const int dim, const int n, Metric& m, NumericVector& _r_k, const double _alpha)
  : DualTreeBoruvka(prune, dim, n, m), alpha(_alpha), r_k(_r_k) {
  N_q_par = N_r_par = NULL;
  if (prune){ bnd_knn = new std::unordered_map<ANNkd_node*, BoundKNN& >(); }
}

// Override the istance calculation to bound the connection radii
ANNdist DTB_CT::computeDistance(const int q_idx, const int r_idx, ANNdist eps1, ANNdist eps2){
  ANNcoord* qq = qtree->pts[q_idx];     // first coord of query point
  ANNcoord* pp = rtree->pts[r_idx];			// first coord of reference point
  ANNdist dist = m_dist(qq, pp);        // Compute (full) metric distance
  return std::max(dist / alpha, std::max(r_k[q_idx], r_k[r_idx])); // Return augmented distance
}