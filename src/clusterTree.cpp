#include <Rcpp.h>
using namespace Rcpp;

// Includes
#include "union_find.h"
#include "utilities.h"

// Computes the connection radius, i.e. the linkage criterion
inline double getConnectionRadius(double r, double radius_i, double radius_j, double alpha, const int type) {

  // Only admit edges with finite weight if allowed
  // double min_radius = std::min(radius_i, radius_j);
  // if (r > min_radius){ return std::numeric_limits<double>::infinity(); }

  switch(type){
  // Robust Single Linkage from 2010 paper
  case 0:
    return std::max(r / alpha, std::max(radius_i, radius_j));
    break;
    // kNN graph from Algorithm 2 from Luxburgs 2014 paper
  case 1:
    // return (radius_i <= (r)) && (radius_j <= (r)) && (dist_ij) <= alpha * std::max(radius_i, radius_j);
    // double max_radius = std::max(radius_i, radius_j);
    return std::max(radius_i, radius_j);
    break;
    // mutual kNN graph from Algorithm 2 from Luxburgs 2014 paper
  case 2:
    return std::min(r / alpha, std::max(radius_i, radius_j));
    break;
  default:
    Rcpp::stop("Not a valid neighborhood query type");
  }
  return false;
}

// Recursively visit the merge matrix to extract an hclust sufficient ordering
void visit(const IntegerMatrix& merge, IntegerVector& order, int i, int j, int& ind) {
  // base case
  if (merge(i, j) < 0) {
    order.at(ind++) = -merge(i, j);
  }
  else {
    visit(merge, order, merge(i, j) - 1, 0, ind);
    visit(merge, order, merge(i, j) - 1, 1, ind);
  }
}

// Build a visually compatible ordering
IntegerVector extractOrder(IntegerMatrix merge){
  IntegerVector order = IntegerVector(merge.nrow()+1);
  int ind = 0;
  visit(merge, order, merge.nrow() - 1, 0, ind);
  visit(merge, order, merge.nrow() - 1, 1, ind);
  return(order);
}

// Regular Kruskal's MST
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


// Algorithm 2: kNN Graph
// r_k := vector of k - 1 nearest neighbors
// knn_indices := id's of knn correspond to r_k
NumericMatrix kruskalsKNN(const NumericVector r_k,
                          const IntegerVector knn_indices,
                          const int n, const int k, const double alpha){

  Rcpp::Rcout << "Beginning kruskals" << std::endl;
  // Set up resulting MST
  NumericMatrix mst = NumericMatrix(n - 1, 3);

  // Get sorted radii
  NumericVector lambda = Rcpp::clone(r_k).sort(false);

  // Get order of original; use R order function to get consistent ordering
  Function order = Function("order"), duplicated = Function("duplicated");
  IntegerVector r_order = as<IntegerVector>(order(r_k)) - 1;

  // Create disjoint-set data structure to track components
  UnionFind components = UnionFind(n);
  int i = 0, crow = 0;
  for (NumericVector::const_iterator r = lambda.begin(); r != lambda.end(); ++r, ++i) {

    // Retrieve index of x_i and x_j
    int x_i = r_order.at(i); //(int) INDEX_FROM_KNN(i+1, k);
    int x_j = int(knn_indices.at(x_i) - 1);
    if (i < 10 || i % 100 == 0) Rcout << "i: " << i << ", " << x_i << " " << x_j << std::endl;
    if (x_i < n && x_j < n){
      if (components.Find(x_i) != components.Find(x_j)){
        mst(crow, _) = NumericVector::create(x_i, x_j, *r);
        crow++;
        components.Union(x_i, x_j);
        if (crow == n) break;
      }
    }
  }

  return(mst);
}

/* primsMST
 * Compute MST using variant of Prim's based on a priority-first search for a dense graph
 */
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


/*
* Compute MST using variant of Prim's, constrained by the radius of the Balls around each x_i.
* Requires several array-type or indicator variables, namely:
* v_selected := array of indicators of spanning tree membership (-1 implies non-membership)
* c_i := index of current (head) node (relative to r_k)
* t_i := index of node to test against (relative to r_k)
* d_i := index of distance from current node to test node (relative to r)
*/
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
          fringe[t_i].weight = priority;
          fringe[t_i].to = c_i; // t_i indexes the 'from' node
        }

        // An edge 'on the fringe' might be less than any of the current nodes weights
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


/* mstToHclust
 * Given a minimum spanning tree of the columnar form (<from>, <to>, <height>), create a valid hclust object
 * using a disjoint-set structure to track components
 * Notes: expects 0-based mst indices, and that all edge indices are 0 <= i < n
 */
// [[Rcpp::export(name = ".mstToHclust")]]
List mstToHclust(NumericMatrix mst, const int n){

  // Extract merge heights and associated order of such heights
  NumericVector dist = mst.column(2);
  IntegerVector height_order = order_(dist) - 1;

  // Set up components
  UnionFind components = UnionFind(n);
  IntegerVector comp_index = rep(0, n);
  std::vector<int> assigned = std::vector<int>(n, 0); // 0 to indicate node has not been assigned

  // Go through each row of MST, tracking component indices
  IntegerMatrix merge = IntegerMatrix(n - 1, 2);
  for (int i = 0; i < n - 1; ++i) {
    if (i % 1000 == 0) Rcpp::checkUserInterrupt();
    NumericVector crow = mst.row(height_order.at(i));
    int from = crow.at(0), to = crow.at(1);
    int from_comp = components.Find(from), to_comp = components.Find(to);
    components.Union(from, to);

    // CASE 1: Two singletons a merging
    if (assigned.at(from) == 0 && assigned.at(to) == 0) { // Rcout << "1";
      merge.row(i) = IntegerVector::create(-(from + 1), -(to + 1));
      comp_index.at(components.Find(from)) = i;
    }
    // CASE 2: One singleton and one component are merging
    else if (assigned.at(from) == 0 || assigned.at(to) == 0) { // Rcout << "2";
      int leaf = assigned.at(from) == 0 ? from : to;
      int comp = assigned.at(from) == 0 ? to_comp : from_comp;
      merge.row(i) = IntegerVector::create(-(leaf+1), comp_index.at(comp)+1);
      comp_index.at(components.Find(leaf)) = i;
    }
    // CASE 3: Two components
    else {
      merge.row(i) = IntegerVector::create(comp_index.at(from_comp)+1, comp_index.at(to_comp)+1);
      comp_index.at(components.Find(from)) = i;
    }
    assigned.at(from) = assigned.at(to) = 1;
  }

  // Extractor merge order and return
  List res = List::create(
    _["merge"] = merge,
    _["height"] = dist[height_order],
    _["order"] = extractOrder(merge),
    _["labels"] = R_NilValue
    //_["mst"] = mst
  );
  res.attr("class") = "hclust";
  return(res);
}

// [[Rcpp::export]]
List clusterTree(const NumericVector dist_x, const NumericVector r_k, const int k, const double alpha = 1.414213562373095,
                 const int type = 0, IntegerVector knn_indices = IntegerVector()) {
  Rcout << "starting cluster tree" << std::endl;
  std::string message = "clusterTree expects a 'dist' object.";
  if (!dist_x.hasAttribute("class") || as<std::string>(dist_x.attr("class")) != "dist") { stop(message); }
  if (!dist_x.hasAttribute("method")) { stop(message); }
  if (!dist_x.hasAttribute("Size")){ stop(message); }
  if (as<std::string>(dist_x.attr("method")) != "euclidean") { warning("RSL expects euclidean distances.");}

  // Number of data points
  const int n = as<int>(dist_x.attr("Size"));

  // Get ordered radii, use R order function to get consistent ordering
  NumericVector r = Rcpp::clone(dist_x);

  //Rcout << "converting to linear" << std::endl;
  //NumericVector r_k_linear = as<NumericVector>(r_k);
  // Run the MST with set parameters


  // Rcout << "converting to linear" << std::endl;
  List res;
  //NumericMatrix mst;
  switch(type){
    case 0:
      res = mstToHclust(primsRSL(r, r_k, n, alpha, 0), n);
      break;
    case 1:
      res = List::create(kruskalsKNN(r_k, knn_indices, n, k, alpha));
      break;
    case 2:
      break;
    default:
      Rcpp::stop("Invalid type supplied.");
      break;
  }

  // Use MST to create condensed hclust hierarchy
  //List res = mstToHclust(mst, n);
  return (res);
}
