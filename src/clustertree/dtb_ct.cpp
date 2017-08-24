#include "dtb_ct.h" // Dual Tree Borvuka definitions

// use default constructor
DTB_CT::DTB_CT(const NumericMatrix& q_x, Metric& m, NumericMatrix& r_x, List& config) : DualTreeBoruvka(q_x, m, r_x, config)
{
  if (!config.containsElementNamed("alpha")) Rcpp::stop("Invalid parameters passed to Dual Tree Boruvka.");
  alpha = (double) config["alpha"];
}

// Override the istance calculation to bound the connection radii
ANNdist DTB_CT::computeDistance(const int q_idx, const int r_idx, ANNdist eps1, ANNdist eps2){
  ANNcoord* qq = qtree->pts[q_idx];     // first coord of query point
  ANNcoord* pp = rtree->pts[r_idx];			// first coord of reference point
  ANNdist dist = m_dist(qq, pp);        // Compute (full) metric distance
  return std::max(dist / alpha, std::max(r_k[q_idx], r_k[r_idx])); // Return augmented distance
}