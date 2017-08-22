#include <Rcpp.h>
using namespace Rcpp;

#include <clustertree/clustertree.h> // Auxiliary C++ extensions related to the clustertree
#include <clustertree/dtb_ct.h> // Dual Tree Boruvka extensions for RSL


// Use the dual tree boruvka approach to compute the cluster tree
// [[Rcpp::export]]
List dtbRSL(const NumericMatrix& x, const NumericVector& r_k, const double alpha, const int type, SEXP metric_ptr){
  const int d = x.ncol();
  const int n = x.nrow();

  // Copy data over to ANN point array
  ANNkd_tree* kd_treeQ, *kd_treeR;
  ANNpointArray x_ann = matrixToANNpointArray(x);

  // Construct the dual tree KNN instance
  Metric& metric = getMetric(metric_ptr);
  DTB_CT dtb_setup = DTB_CT(true, d, n, metric, r_k, alpha);

  // Construct the tree
  ANNkd_tree* kd_tree = dtb_setup.ConstructTree(x_ann, x.nrow(), x.ncol(), 30, ANN_KD_SUGGEST);

  // With the tree(s) created, setup DTB-specific bounds, assign trees, etc.
  dtb_setup.setup(kd_tree, kd_tree);

  // Debug mode: Print the tree
  UTIL(dtb.PrintTree((ANNbool) true, true))

  // Run the dual tree boruvka algorithm (w/ augmented distance function)
  List mst = dtb_setup.DTB(x);

  return mst;
}

// [[Rcpp::export]]
List clusterTree(const NumericVector dist_x, const NumericVector r_k, const int k, const double alpha = 1.414213562373095,
                 const int type = 0) {
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
  List res = mstToHclust(mst);
  return (res);
}