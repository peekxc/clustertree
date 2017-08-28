#include <Rcpp.h>
using namespace Rcpp;

#include "hclust_util.h" // Hclust extensions
#include "ANN_util.h" // matrixToANNpointArray
#include "metrics.h"

// Computes the connection radius, i.e. the linkage criterion
double getConnectionRadius(double dist_ij, double radius_i, double radius_j, double alpha, const int type) {

  // Only admit edges with finite weight if the neighborhood radii allow
  // Note that RSL will always form a complete hierarchy, so returning the numerical
  // limits maximum isn't necessary.
  double R;
  switch(type){
  // Robust Single Linkage from 2010 paper
  case 0:
    return std::max(dist_ij / alpha, std::max(radius_i, radius_j));
    break;
    // kNN graph from Algorithm 2 from Luxburgs 2014 paper
  case 1:
    R = alpha * std::max(radius_i, radius_j);
    return dist_ij <= R ? R : std::numeric_limits<double>::max();
    break;
    // mutual kNN graph from Algorithm 2 from Luxburgs 2014 paper
  case 2:
    R = alpha * std::min(radius_i, radius_j);
    return dist_ij <= R ? R : std::numeric_limits<double>::max();
    break;
  default:
    Rcpp::stop("Not a valid neighborhood query type");
  }
  return std::numeric_limits<double>::max();
}

// // Use the dual tree boruvka approach to compute the cluster tree
// List dtbRSL(const NumericMatrix& x, const NumericVector& r_k, const double alpha, const int type, SEXP metric_ptr){
//   const int d = x.ncol();
//   const int n = x.nrow();
//
//   // Copy data over to ANN point array
//   // ANNkd_tree* kd_treeQ, *kd_treeR;
//   ANNpointArray x_ann = matrixToANNpointArray(x);
//
//   // Construct the dual tree KNN instance
//   Metric& metric = getMetric(metric_ptr);
//   const NumericMatrix& q_x, Metric& m, NumericMatrix& r_x = emptyMatrix, List config = List::create()
//   DTB_CT dtb_setup = DTB_CT(x, metric, emptyMatrix, alpha);
//
//   // Construct the tree
//   ANNkd_tree* kd_tree = dtb_setup.ConstructTree(x_ann, x.nrow(), x.ncol(), 30, ANN_KD_SUGGEST);
//
//   // With the tree(s) created, setup DTB-specific bounds, assign trees, etc.
//   dtb_setup.setup(kd_tree, kd_tree);
//
//   // Run the dual tree boruvka algorithm (w/ augmented distance function)
//   List mst = dtb_setup.DTB(x);
//
//   return mst;
// }

/*
 * Compute MST using variant of Prim's, constrained by the radius of the Balls around each x_i.
 * Requires several array-type or indicator variables, namely:
 * v_selected := array of indicators of spanning tree membership (-1 implies non-membership)
 * c_i := index of current (head) node (relative to r_k)
 * t_i := index of node to test against (relative to r_k)
 * d_i := index of distance from current node to test node (relative to r)
 */
// [[Rcpp::export]]
NumericMatrix primsRSL(const NumericVector r, const NumericVector r_k, const int n, const double alpha, const int type){
  // Set up resulting MST
  NumericMatrix mst = NumericMatrix(n - 1, 3);

  // Data structures for prims
  std::vector<int> v_selected = std::vector<int>(n, -1); // -1 to indicate node is not in MST
  std::vector<edge> fringe = std::vector<edge>(n, edge(-1, std::numeric_limits<double>::infinity()));

  double min = std::numeric_limits<double>::infinity(), priority = 0.0;
  int c_i = 0, min_id = n - 1;
  for (int n_edges = 0; n_edges < n - 1; n_edges++) {
    if (n_edges % 1000 == 0) Rcpp::checkUserInterrupt();
    min = std::numeric_limits<double>::infinity(); // Reset so new edge is always chosen

    // Compare all the new edge weights w/ the "current best" edge weights
    for (int t_i = 0; t_i < n; ++t_i) {
      // Graph subset step: ensure both c_i and t_i have neighborhood radii at least as big as the current radius
      if (t_i == c_i) continue;
      int d_i = t_i > c_i ? INDEX_TF(n, c_i, t_i) : INDEX_TF(n, t_i, c_i); // bigger index always on the right

      // MST step, make sure node isn't already in the spanning tree
      if (v_selected[t_i] < 0) {

        // Generic RSL step
        priority = getConnectionRadius(r[d_i], r_k[c_i], r_k[t_i], alpha, type);
        if (priority < fringe[t_i].weight) { // using max above implicitly ensures c_i and t_i are connected
          // Rcout << "Updating fringe " << t_i << "w/ radii (" << r_k[c_i] << ", " << r_k[t_i] << ")" << std::endl;
          // Rcout << "current: " << c_i << ", to: " << t_i << std::endl;
          fringe[t_i].weight = priority;
          fringe[t_i].to = c_i; // t_i indexes the 'from' node
        }

        // An edge 'on the fringe' might be less than any of the current nodes weights
        if (fringe[t_i].weight < min) {
          min = fringe[t_i].weight, min_id = t_i;
        }
      }
    }
    // Rcout << "Adding edge: (" << min_id << ", " << c_i << ") [" << min << "]" << std::endl;
    mst(n_edges, _) = NumericVector::create(min_id, c_i, min);
    v_selected[c_i] = 1;
    c_i = min_id;
  }
  return(mst);
}

// [[Rcpp::export]]
NumericMatrix naive_clustertree(const NumericVector& dist_x, const NumericVector& r_k, const int k, const double alpha, const int type = 0) {
  std::string message = "naive_clustertree expects a 'dist' object.";
  if (!dist_x.hasAttribute("class") || as<std::string>(dist_x.attr("class")) != "dist") { stop(message); }
  if (!dist_x.hasAttribute("method")) { stop(message); }
  if (!dist_x.hasAttribute("Size")){ stop(message); }

  // Get sorted radii
  const int n = as<int>(dist_x.attr("Size"));
  NumericVector sorted_x = Rcpp::clone(dist_x).sort(false);
  NumericMatrix mst = NumericMatrix(n - 1, 3);

  // Get order of original; use R order function to get consistent ordering
  Function order = Function("order");
  IntegerVector r_order = as<IntegerVector>(order(dist_x)) - 1;

  // Create disjoint-set data structure to track components
  UnionFind components = UnionFind(n);
  int i = 0, crow = 0;
  double r = 0;
  for (NumericVector::const_iterator dist_ij = sorted_x.begin(); dist_ij != sorted_x.end(); ++dist_ij, ++i) {
    // Retrieve index of x_i and x_j
    int x_i = INDEX_TO(i, n), x_j = INDEX_FROM(i, n, x_i);
    r = *dist_ij / alpha;
    Rcout << "Comparing: " << x_i << ", " << x_j << " (r = " << r << ")" << std::endl;
    if (r_k[x_i] <= r && r_k[x_j] <= r){ // if admitted
      Rcout << "Checking: " << x_i << "(r_k = " << r_k[x_i] << "), " << x_j << "(r_k = " << r_k[x_j] << ") against " << r << std::endl;
      if (components.Find(x_i) != components.Find(x_j)){
        Rcout << "Connecting: " << x_i << ", " << x_j << " @r = " << r << std::endl;
        mst(crow++, _) = NumericVector::create(x_i, x_j, r); // Use r to index cluster tree
      }
      components.Union(x_i, x_j);
    }
  }
  return mst;
}
