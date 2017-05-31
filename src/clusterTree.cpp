#include <Rcpp.h>
using namespace Rcpp;

// Includes
#include "union_find.h"
#include "utilities.h"

// Allows indexing lower triangular (dist objects)
#include <math.h>
#define INDEX_TF(N,to,from) (N)*(to) - (to)*(to+1)/2 + (from) - (to) - (1)
#define INDEX_TO(k, n) n - 2 - floor(sqrt(-8*k + 4*n*(n-1)-7)/2.0 - 0.5)
#define INDEX_FROM(k, n, i) k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2

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
List checkKruskals(const NumericVector x, const NumericVector r_k, const int k, const double alpha = 1.414213562373095){

  // Number of data points
  const int n = as<int>(x.attr("Size"));

  // Get sorted radii
  NumericVector radii = Rcpp::clone(x).sort(false);// * alpha;

  // Get order of original; use R order function to get consistent ordering
  Function order = Function("order");
  IntegerVector r_order = as<IntegerVector>(order(x)) - 1;

  // Set up connected components (initialized as singletons)
  UnionFind components = UnionFind(n);

  // Set up vector whose values, indexed by component, represent last merge step the component was a part of
  IntegerVector component_index = IntegerVector(n);

  // Set up vector representing inclusion and admittance criteria
  std::vector<bool> admitted = std::vector<bool>(n, false);
  std::vector<bool> included = std::vector<bool>(n, false);
  //IntegerVector included = Rcpp::seq(0, n - 1);

  // Admission tracker
  NumericVector sorted_r_k = clone(r_k).sort(false);
  IntegerVector order_r_k = as<IntegerVector>(order(r_k)) - 1;
  int c_rk = 0;

  // Iterate r from 0 \to \infty
  IntegerMatrix order_ft = IntegerMatrix(radii.length(), 2);
  int i = 0, c = 0;
  IntegerMatrix merge = IntegerMatrix(n - 1, 2);
  IntegerMatrix merge_original = IntegerMatrix(n - 1, 2);
  NumericVector height = NumericVector(n - 1);
  for(NumericVector::iterator r = radii.begin(); r != radii.end(); ++r, ++i){
    int to = INDEX_TO(r_order.at(i), n), from = INDEX_FROM(r_order.at(i), n, to);
    int from_comp = components.Find(from), to_comp = components.Find(to);

    // Counter to current lowest neighborhood radius; progressively sets admission of points
    while(c_rk < n && *r > sorted_r_k.at(c_rk)) admitted.at(int(order_r_k.at(c_rk++))) = true;

    // Construct a graph G_r with nodes \{ x_i : r_k(x_i) \leq r \}.
    if (admitted.at(from) && admitted.at(to)){
      // Include edge (x_i, x_j) if \lVert x_i - x_j \rVert \leq \alpha r.
      if (components.Find(from) != components.Find(to)){
        components.Union(from, to);
        height.at(c) = (*r) * alpha;
        merge_original(c, _) = IntegerVector::create(from, to);

        // Hclust class requirements: when building the merge matrix, positive merge indices represent agglomerations,
        // and negative indicate singletons. This requires more indexing vectors to enable tracking component
        // membership changes
        if (!included.at(from) && !included.at(to)){
          Rcout << "Merging 0: " << -from << " and " << -to << std::endl;
          merge(c++, _) = IntegerVector::create(-(from + 1), -(to + 1));
          included.at(from) = included.at(to) = true;
          //component_index.at(components.Find(from)) = c;
        } else if (!included.at(from) || !included.at(to)){
          int leaf = (!included.at(from)) ? from : to;
          int comp = (!included.at(from)) ? to_comp : from_comp;
          Rcout << "Merging 1: " << -leaf << " and " << comp << std::endl;
          merge(c++, _) = IntegerVector::create(-(leaf + 1), component_index.at(comp));
          included.at(leaf) = true;
          //component_index.at(components.Find(leaf)) = c;
        } else {
          Rcout << "Merging 2: " << from_comp << " and " << to_comp << std::endl;
          merge(c++, _) = IntegerVector::create(component_index.at(from_comp), component_index.at(to_comp));
          //component_index.at(components.Find(from)) = c;
        }
        Rcout << "updating component index: " << from_comp << ", " << to_comp << std::endl;
        component_index.at(from_comp) = component_index.at(to_comp) = c;
        // if (c == 2) break;
        //merge(c++, _) = IntegerVector::create(from, to);
      }
    }
    order_ft(i, _) = IntegerVector::create(from, to);
  }
  Rcout << "here" << std::endl;
  IntegerVector final_component = IntegerVector(n, 0);
  for (int i = 0; i < n; ++i){
    final_component.at(i) = components.Find(i);
  }


  // Extractor merge order and return
  List res = List::create(
    _["merge"] = merge,
    _["height"] = height,
    _["order"] = extractOrder(merge),
    _["labels"] = R_NilValue
  );
  res.attr("class") = "hclust";

  return(List::create(_["merge"] = merge, _["merge_original"] = merge_original, _["height"] = height, _["ind"] = order_ft, _["fc"] = final_component,
                      _["admitted"] = wrap(admitted), _["comp_index"] = component_index, _["hc"] = res));
}

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
