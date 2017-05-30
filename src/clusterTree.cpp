#include <Rcpp.h>
using namespace Rcpp;

// Includes
#include "union_find.h"
#include "utilities.h"

// Allows indexing lower triangular (dist objects)
#define INDEX_TF(N,to,from) (N)*(to) - (to)*(to+1)/2 + (from) - (to) - (1)

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

// Structures to do priority queue
struct edge
{
  unsigned int to;
  double weight;
  edge(int to_id, double cost) : to(to_id), weight(cost) { }
};
struct compare_edge
{
  bool operator()(const edge& e1, const edge& e2) const
  { return e1.weight > e2.weight; }
};

// [[Rcpp::export]]
List clusterTree(const NumericVector x, const NumericVector r_k, const int k, const double alpha = 1.414213562373095) {
  std::string message = "clusterTree expects a 'dist' object.";
  if (!x.hasAttribute("class") || as<std::string>(x.attr("class")) != "dist") { stop(message); }
  if (!x.hasAttribute("method")) { stop(message); }
  if (!x.hasAttribute("Size")){ stop(message); }
  if (as<std::string>(x.attr("method")) != "euclidean") { warning("RSL expects euclidean distances.");}

  // Number of data points
  const int n = as<int>(x.attr("Size"));

  // Get ordered radii, use R order function to get consistent ordering
  NumericVector r = Rcpp::clone(x);// * alpha;

  // Resulting MST
  NumericMatrix mst = NumericMatrix(n - 1, 3);

  // Data structures for prims
  std::vector<int> v_selected = std::vector<int>(n, -1); // -1 to indicate node is not in MST
  std::vector<edge> fringe = std::vector<edge>(n, edge(-1, std::numeric_limits<double>::infinity()));

  /* PHASE 1
   * Compute MST using variant of Prim's, constrained by the radius of the Balls around each x_i.
   * Requires several array-type or indicator variables, namely:
   * v_selected := array of indicators of spanning tree membership (-1 implies non-membership)
   * c_i := index of current (head) node (relative to r_k)
   * t_i := index of node to test against (relative to r_k)
   * d_i := index of distance from current node to test node (relative to r)
   */
  double min = std::numeric_limits<double>::infinity(), priority = 0.0;
  int c_i = 0, min_id = n - 1;
  for (int n_edges = 0; n_edges < n - 1; n_edges++) {
    if (n_edges % 1000 == 0) Rcpp::checkUserInterrupt();
    min = std::numeric_limits<double>::infinity(); // Reset so new edge is always chosen
    // Compare all the new edge weight w/ the "current best" edge weights
    for (int t_i = 0; t_i < n; ++t_i) {
      // Graph subset step: ensure both c_i and t_i have neighborhood radii at least as big as the current radius
      if (t_i == c_i) continue;
      int d_i = t_i > c_i ? INDEX_TF(n, c_i, t_i) : INDEX_TF(n, t_i, c_i); // bigger index always on the right
      //if (r[d_i] < std::max(r_k[c_i], r_k[t_i])) continue;

      // MST step, make sure node isn't already in the spanning tree
      if (v_selected[t_i] < 0) {

        // RSL step: connect c_i and t_i if ||c_i - y_i|| <= \alpha * r
        priority = std::max(r[d_i] * alpha, std::max(r_k[c_i], r_k[t_i]));
        if (priority < fringe[t_i].weight) { // using max above implicitly ensures c_i and t_i are connected
          fringe[t_i].weight = priority;
          fringe[t_i].to = c_i; // t_i indexes the 'from' node
        }

        if (fringe[t_i].weight < min) { // an edge 'on the fringe' might be less than any of the current nodes weights
          min = fringe[t_i].weight, min_id = t_i;
        }
      }
    }

    mst(n_edges, _) = NumericVector::create(min_id, c_i, min);
    v_selected[c_i] = 1;
    c_i = min_id;
  }
  return(List::create(mst, v_selected));

  /* PHASE 2
   * Sort MST and use disjoint-set to track components
   */
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

    // Rcout << "(" << i << ") from: " << from << ", to: " << to << std::endl;
    // Rcout << "(" << i << ") from_comp: " << from_comp << ", to_comp: " << to_comp << std::endl;
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
    else { // Rcout << "3";
      merge.row(i) = IntegerVector::create(comp_index.at(from_comp)+1, comp_index.at(to_comp)+1);
      comp_index.at(components.Find(from)) = i;
    }
    assigned.at(from) = assigned.at(to) = 1;
    // if (i == 1) return(List::create(assigned, comp_index, merge, mst, height_order));
  }

  // Extractor merge order and return
  List res = List::create(
    _["merge"] = merge,
    _["height"] = dist[height_order],
    _["order"] = extractOrder(merge),
    _["labels"] = R_NilValue,
    _["mst"] = mst
  );
  res.attr("class") = "hclust";
  return (res);
}
