#include "dtb_ct.h" // Dual Tree Borvuka definitions

// use default constructor
DTB_CT::DTB_CT(const NumericMatrix& q_x, Metric& m, NumericMatrix& r_x, List& config) : DualTreeBoruvka(q_x, m, r_x, config)
{
  if (!config.containsElementNamed("alpha")) Rcpp::stop("Invalid parameters passed to Cluster tree MST.");
  if (!config.containsElementNamed("R")) Rcpp::stop("Invalid parameters passed to Cluster tree MST.");
  if (!config.containsElementNamed("estimator")) Rcpp::stop("Invalid parameters passed to Cluster tree MST.");
  alpha = (double) config["alpha"];
  r_k = as<NumericVector>(config["R"]);
  type = config["estimator"]; // 0 for RSL, 1 for knn, 2 for Mutual KNN
}

// Override the istance calculation to bound the connection radii
ANNdist DTB_CT::computeDistance(const int q_idx, const int r_idx, ANNdist eps1, ANNdist eps2){
  //R_INFO("Computing radii-constrained distance\n");
  ANNcoord* qq = qtree->pts[q_idx];     // first coord of query point
  ANNcoord* pp = rtree->pts[r_idx];			// first coord of reference point
  ANNdist dist = m_dist(qq, pp);        // Compute (full) metric distance
  return getConnectionRadius(dist, r_k[q_idx], r_k[r_idx], alpha, type); // Return augmented distance
}


// Cluster tree connection radius function.
// To allow the dual tree MST algorithm to finish, the maximum return value is the maximum allowed by the precision of
// ANNdist. Inf will cause the program to never cease, and the algorithm will keep attempting to find the minimum
// spanning edge, without ever finding one.
ANNdist DTB_CT::getConnectionRadius(double dist_ij, double radius_i, double radius_j, double alpha, const int type){

  // Only admit edges with finite weight if the neighborhood radii allow
  // Note that RSL will always form a complete hierarchy, so returning the numerical
  // limits maximum isn't necessary.
  ANNdist R;
  switch(type){
  // Robust Single Linkage from 2010 paper
  case 0:
    return std::max(dist_ij / alpha, std::max(radius_i, radius_j));
    break;
    // kNN graph from Algorithm 2 from Luxburgs 2014 paper
  case 1:
    R = alpha * std::max(radius_i, radius_j);
    return dist_ij <= R ? R : std::numeric_limits<ANNdist>::max();
    break;
    // mutual kNN graph from Algorithm 2 from Luxburgs 2014 paper
  case 2:
    R = alpha * std::min(radius_i, radius_j);
    return dist_ij <= R ? R : std::numeric_limits<ANNdist>::max();
    break;
  default:
    Rcpp::stop("Not a valid neighborhood query type");
  }
  return std::numeric_limits<ANNdist>::max();
}