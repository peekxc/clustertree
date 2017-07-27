#include <Rcpp.h>
using namespace Rcpp;

#include <clustertree/clustertree.h>

// [[Rcpp::export]]
List clusterTree(const NumericVector dist_x, const NumericVector r_k, const int k, const double alpha = 1.414213562373095,
                 const int type = 0, IntegerVector knn_indices = IntegerVector()) {
  std::string message = "clusterTree expects a 'dist' object.";
  if (!dist_x.hasAttribute("class") || as<std::string>(dist_x.attr("class")) != "dist") { stop(message); }
  if (!dist_x.hasAttribute("method")) { stop(message); }
  if (!dist_x.hasAttribute("Size")){ stop(message); }
  if (as<std::string>(dist_x.attr("method")) != "euclidean") { warning("RSL expects euclidean distances.");}

  // Number of data points
  const int n = as<int>(dist_x.attr("Size"));

  // Get ordered radii, use R order function to get consistent ordering
  NumericVector r = Rcpp::clone(dist_x);

  // Run the MST with set parameters
  NumericMatrix mst = primsRSL(r, r_k, n, alpha, type);

  // Convert to HCLUST object
  List res = mstToHclust(mst, n);
  return (res);
}