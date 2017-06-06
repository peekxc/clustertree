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

// Retrieve the full integer vector representing component membership (by reference)
void CC(UnionFind& uf, IntegerVector& cc, std::vector<bool>& admitted){
  for (int i = 0; i < uf.size; ++i){
    cc[i] = admitted[i] ? uf.Find(i) : 0;
  }
}

// Retrieve the full integer vector representing component membership (by reference)
void CC(UnionFind& uf, IntegerVector& cc, std::vector<int>& admitted){
  for (std::vector<int>::iterator it = admitted.begin(); it < admitted.end(); ++it){
    cc[*it] = uf.Find(*it);
  }
}

// Retrieve the full integer vector representing component membership
IntegerVector CC(UnionFind& uf){
  IntegerVector cc = Rcpp::no_init(uf.size);
  for (int i = 0; i < uf.size; ++i){ cc[i] = uf.Find(i); }
  return(cc);
}

void printCC(UnionFind& uf, std::vector<bool>& admitted){
  // for (std::vector<int>::iterator it = admitted.begin(); it < admitted.end(); ++it){
  //   cc[*it] = uf.Find(*it);
  // }
  for (int i = 0; i < uf.size; ++i){
    Rcout << (admitted[i] ? uf.Find(i) : 0) << " ";
  }
}

/* -- CC_changed --
 * This function controls for the detection of a change between two sets of connected
 * components, restricted to only look for changes between points that have been admitted into
 * G_r.
 */
bool CC_changed(IntegerVector current_CC, IntegerVector previous_CC, std::vector<int>& admitted){
  if (admitted.size() < 2) return(false);
  for (std::vector<int>::iterator it = admitted.begin(); it != admitted.end(); ++it){
    if (current_CC[*it] != previous_CC[*it]) {// && contains(admitted, current_CC[*it])){
      return(true);
    }
  }
  return(false);
}


// [[Rcpp::export]]
List checkKruskals(const NumericVector x, const NumericVector r_k, const int k, const double alpha = 1.414213562373095){

  // Number of data points
  const int n = as<int>(x.attr("Size"));

  // Get sorted radii
  NumericVector radii = Rcpp::clone(x).sort(false) * alpha;
  NumericVector l2 = Rcpp::clone(x).sort(false) ;

  // Get order of original; use R order function to get consistent ordering
  Function order = Function("order");
  IntegerVector r_order = as<IntegerVector>(order(x)) - 1;

  // Set up connected components (initialized as singletons)
  UnionFind components = UnionFind(n);

  // Set up vector representing inclusion and admittance criteria
  // admittance := only consider x_i s.t. \{ x_i : r_k(x_i) \leq r \}
  // included := 0 if x_i is a singleton, 1 if x_i is connected to some component
  std::vector<bool> admitted = std::vector<bool>(n, false);
  std::vector<int> admitted2 = std::vector<int>();
  std::vector<bool> included = std::vector<bool>(n, false);

  // Admission tracker; rather than check admissability of every point at every radii,
  // sort by increasing radius and make one comparison each loop
  int c_rk = 0;
  NumericVector sorted_r_k = clone(r_k).sort(false);
  IntegerVector order_r_k = as<IntegerVector>(order(r_k)) - 1;


  IntegerVector index_conjoined = IntegerVector(n, 0);

  int i = 0, c = 0;
  // Data structures needed to build the 'hclust' object
  IntegerVector component_index = IntegerVector(n);

  // Set up vector whose values, indexed by component, represent last merge step the component was a part of
  IntegerMatrix merge = IntegerMatrix(n - 1, 2);
  NumericVector height = NumericVector(n - 1);
  IntegerMatrix merge_order = IntegerMatrix(radii.length(), 2);

  List cc_list = List(n, IntegerVector());

  IntegerVector previous_CC = seq(0, n - 1);
  IntegerVector current_CC = clone(previous_CC);

  NumericVector R = NumericVector(n);

  // Iterate r from 0 \to \infty
  for(NumericVector::iterator r = radii.begin(), dist_ij = l2.begin(); r != radii.end(); ++r, ++i, ++dist_ij){
    int to = INDEX_TO(r_order.at(i), n), from = INDEX_FROM(r_order.at(i), n, to);
    int from_comp = components.Find(from), to_comp = components.Find(to);
    merge_order(i, _) = IntegerVector::create(from, to);

    // Sets admission status of points
    while(c_rk < n && sorted_r_k.at(c_rk) <= (*dist_ij)) {
      admitted.at(int(order_r_k.at(c_rk))) = true;
      admitted2.push_back(int(order_r_k.at(c_rk)));
      //Rcout << i << ": Set " << int(order_r_k.at(c_rk)) << " admitted == TRUE (@radius == " << *dist_ij << ")" << std::endl;
      //Rcout << "HEAD: " << sorted_r_k.at(c_rk) << ", " << sorted_r_k.at(c_rk+1) << ", " << sorted_r_k.at(c_rk+2) << std::endl;
      c_rk++;
      if (i == 33 || i == 34 || i == 35){
        Rcout << "N admitted: " << admitted2.size() << std::endl;
      }
    }

    // Union the points, (potentially) mutating the current set of connected components
    //if (admitted.at(from) && admitted.at(to)){
      components.Union(from, to);
    //}



    // Update current connected components
    //CC(components, current_CC, admitted);
    CC(components, current_CC, admitted2);

    if (i == 28 || i == 29 || i == 30){
      Rcout << i << ": ";
      printCC(components, admitted);
      Rcout << std::endl;
      // Rcout << "from: " << from << ", to: " << to << std::endl;
      // Rcout << "Admitted: ";
      // for (int ii = 0; ii < n; ++ii){
      //   Rcout << admitted.at(ii);
      // }
      // Rcout << std::endl << "radius: " << *r << std::endl;
      // Rcout << "From comp: " << from_comp << ", To comp: " << to_comp << std::endl;
      // Rcout << "Admitted? " << admitted.at(from) << " and " << admitted.at(to) << std::endl;
      Rcout << "Component change? " << CC_changed(current_CC, previous_CC, admitted2) << std::endl;
    }

    // Construct a graph G_r with nodes \{ x_i : r_k(x_i) \leq r \}.
    // if (admitted.at(from) && admitted.at(to)){

      // Include edge (x_i, x_j) if \lVert x_i - x_j \rVert \leq \alpha r.
      // If the edge inclusion results in a component change, record the change
      if (CC_changed(current_CC, previous_CC, admitted2)){
        Rcout << i << ": Including edge: " << from << " -- " << to << " | ";
        printCC(components, admitted);
        Rcout << std::endl;
        index_conjoined.at(c) = i;
        height.at(c) = (*r);
        R.at(c) = (*dist_ij);

        // Save current connected components
        cc_list.at(c) = clone(current_CC);

        // Hclust class requirements: when building the merge matrix, positive merge indices represent agglomerations,
        // and negative indicate singletons. This requires more indexing vectors to enable tracking component
        // membership changes
        if (!included.at(from) && !included.at(to)) {
          merge(c++, _) = IntegerVector::create(-(from + 1), -(to + 1));
          included.at(from) = included.at(to) = true;
        } else if (!included.at(from) || !included.at(to)){
          int leaf = (!included.at(from)) ? from : to;
          int comp = (!included.at(from)) ? to_comp : from_comp;
          merge(c++, _) = IntegerVector::create(-(leaf + 1), component_index.at(comp));
          included.at(leaf) = true;
        } else {
          merge(c++, _) = IntegerVector::create(component_index.at(from_comp), component_index.at(to_comp));
        }
        component_index.at(from_comp) = component_index.at(to_comp) = c;
      }
      previous_CC = clone(current_CC);
      // components.Union(from, to);
    }
  //}
  // IntegerVector final_component = IntegerVector(n, 0);
  // for (int i = 0; i < n; ++i){
  //   final_component.at(i) = components.Find(i);
  // }


  // Extractor merge order and return
  List res = List::create(
    //_["merge"] = merge,
    //_["height"] = height,
    //_["order"] = extractOrder(merge),
    _["labels"] = R_NilValue,
    _["i"] = index_conjoined,
    _["mo"] = merge_order,
    _["cc"] = cc_list,
    _["dist_ij"] = R
  );
  res.attr("class") = "hclust";

  return(res);
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
  //return(List::create(mst, v_selected));

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
