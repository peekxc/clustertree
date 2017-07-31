#include <Rcpp.h>
using namespace Rcpp;

#include "utilities.h" // Indexing macros
#include "ANN/ANN_util.h" // matrixToANNpointArray
#include "DT/dt.h"
#include "DT/MST/dtb.h"
#include "DT/structures/union_find.h" // disjoint set data structure

/* kruskalsMST
 * Compute MST using Kruskals algorithm w/ disjoint set */
// [[Rcpp::export]]
NumericMatrix kruskalsMST(const NumericVector dist_x){
  std::string message = "kruskalsMST expects a 'dist' object.";
  if (!dist_x.hasAttribute("class") || as<std::string>(dist_x.attr("class")) != "dist") { stop(message); }
  if (!dist_x.hasAttribute("method")) { stop(message); }
  if (!dist_x.hasAttribute("Size")){ stop(message); }

  // Set up resulting MST
  const int n = as<int>(dist_x.attr("Size"));
  NumericMatrix mst = NumericMatrix(n - 1, 3);

  // Get sorted radii
  NumericVector sorted_x = Rcpp::clone(dist_x).sort(false);

  // Get order of original; use R order function to get consistent ordering
  Function order = Function("order");
  IntegerVector r_order = as<IntegerVector>(order(dist_x)) - 1;

  // Create disjoint-set data structure to track components
  UnionFind components = UnionFind(n);
  int i = 0, crow = 0;
  for (NumericVector::const_iterator dist_ij = sorted_x.begin(); dist_ij != sorted_x.end(); ++dist_ij, ++i) {
    // Retrieve index of x_i and x_j
    int x_i = INDEX_TO(i, n), x_j = INDEX_FROM(i, n, x_i);
    if (components.Find(x_i) != components.Find(x_j)){
      mst(crow++, _) = NumericVector::create(x_i, x_j, *dist_ij);
      components.Union(x_i, x_j);
    }
  }
  return(mst);
}

/* primsMST
 * Compute MST using variant of Prim's based on a priority-first search for a dense graph */
// [[Rcpp::export]]
NumericMatrix primsMST(const NumericVector dist_x){
  std::string message = "primsMST expects a 'dist' object.";
  if (!dist_x.hasAttribute("class") || as<std::string>(dist_x.attr("class")) != "dist") { stop(message); }
  if (!dist_x.hasAttribute("method")) { stop(message); }
  if (!dist_x.hasAttribute("Size")){ stop(message); }

  // Set up resulting MST
  const int n = as<int>(dist_x.attr("Size"));
  NumericMatrix mst = NumericMatrix(n - 1, 3);

  // Data structures for prims
  std::vector<int> v_selected = std::vector<int>(n, -1); // -1 to indicate node is not in MST
  std::vector<edge> fringe = std::vector<edge>(n, edge(-1, std::numeric_limits<double>::infinity()));

  double min, cedge_weight = 0.0;
  int c_i = 0, min_id = n - 1;
  for (int n_edges = 0; n_edges < n - 1; n_edges++) {
    if (n_edges % 1000 == 0) Rcpp::checkUserInterrupt();
    min = std::numeric_limits<double>::infinity(); // Reset so new edge is always chosen

    // Compare all the new edge weights w/ the "current best" edge weights
    for (int t_i = 0; t_i < n; ++t_i) {
      if (t_i == c_i) continue;
      int d_i = t_i > c_i ? INDEX_TF(n, c_i, t_i) : INDEX_TF(n, t_i, c_i); // bigger index always on the right

      // MST step, make sure node isn't already in the spanning tree
      if (v_selected[t_i] < 0) {

        // Updates edges on the fringe with lower weights if detected
        cedge_weight = dist_x[d_i];
        if (cedge_weight < fringe[t_i].weight) {
          fringe[t_i].weight = cedge_weight;
          fringe[t_i].to = c_i; // t_i indexes the 'from' node
        }

        // If edge 'on the fringe' is less than any of the current (head) nodes weight,
        // change the head-facing node to the edge of the fringe
        if (fringe[t_i].weight < min) {
          min = fringe[t_i].weight, min_id = t_i;
        }
      }
    }

    mst(n_edges, _) = NumericVector::create(min_id, c_i, min);
    v_selected[c_i] = 1;
    c_i = min_id;
  }
  return(mst);
}

// [[Rcpp::export]]
List dtb(NumericMatrix x, const int bkt_size = 30, bool prune = true) {
  // Copy data over to ANN point array
  ANNkd_tree* kd_treeQ, *kd_treeR;
  ANNpointArray x_ann = matrixToANNpointArray(x);

  // Construct the dual tree KNN instance
  DualTreeBoruvka dtb = DualTreeBoruvka(prune, x.ncol(), x.nrow());

  // Construct the tree
  ANNkd_tree* kd_tree = dtb.ConstructTree(x_ann, x.nrow(), x.ncol(), bkt_size, ANN_KD_SUGGEST);

  // With the tree(s) created, setup DTB-specific bounds
  dtb.setup(kd_tree, kd_tree);
}
