// Includes
#include "clustertree.h"

// Computes the connection radius, i.e. the linkage criterion
inline double getConnectionRadius(double dist_ij, double radius_i, double radius_j, double alpha, const int type) {

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

// Algorithm 2: kNN Graph
// r_k := vector of k - 1 nearest neighbors
// knn_indices := id's of knn correspond to r_k
List kruskalsKNN(const NumericVector dist_x,
                          const NumericVector r_k,
                          const IntegerVector knn_indices,
                          const int n, const int k, const double alpha){

  Rcpp::Rcout << "Beginning kruskals" << std::endl;
  // Set up resulting MST
  NumericMatrix mst = NumericMatrix(n - 1, 3);

  // Get sorted radii
  NumericVector lambda = Rcpp::clone(r_k).sort(false);

  // Get order of original; use R order function to get consistent ordering
  Function order = Function("order"), duplicated = Function("duplicated");
  IntegerVector rk_order = as<IntegerVector>(order(r_k)) - 1;
  LogicalVector admitted = LogicalVector(n, false);

  IntegerVector x_order = as<IntegerVector>(order(dist_x)) - 1;
  NumericVector inc_dist = Rcpp::clone(dist_x).sort(false);

  // Create disjoint-set data structure to track components
  UnionFind components = UnionFind(n);
  UnionFind prev_components = UnionFind(n);
  UnionFind dist_components = UnionFind(n);

  bool connect_next = false;
  int i = 0, cr_k = 0, crow = 0, px_i = INDEX_TO(x_order.at(0), n), px_j = INDEX_FROM(x_order.at(0), n, px_i);
  for(NumericVector::iterator dist_ij = inc_dist.begin(); dist_ij != inc_dist.end(); ++dist_ij, ++i){
    int x_j = INDEX_TO(x_order.at(i), n), x_i = INDEX_FROM(x_order.at(i), n, x_j);
    if (dist_components.Find(x_i) != dist_components.Find(x_j)) dist_components.Union(x_i, x_j);
    while(cr_k < n - 1 && r_k.at(rk_order.at(cr_k)) <= (*dist_ij)) { admitted.at(int(rk_order.at(cr_k++))) = true; }
    if (connect_next){
      prev_components.Union(px_i, px_j);
      connect_next = false;
    }
    if ((*dist_ij)/alpha <= lambda.at(cr_k) && admitted.at(x_i) && admitted.at(x_j)){
        //IntegerVector CCs = components.getCC();
        components.merge(dist_components, admitted);
        if (components != prev_components)
        {
          mst(crow, _) = NumericVector::create(x_i, x_j, *dist_ij);
          crow++;
          components.Union(x_i, x_j);
          px_i = x_i, px_j = x_j;
          connect_next = true;
          if (crow == n) break;
        }
    }
  }

  // // Create disjoint-set data structure to track components
  // UnionFind components = UnionFind(n);
  // i = 0;
  // int crow = 0;
  // for (NumericVector::const_iterator r = lambda.begin(); r != lambda.end(); ++r, ++i) {
  //   if (i % 100 == 0) Rcpp::checkUserInterrupt();
  //
  //   // Point x_i becomes admitted into the graph G_r
  //   admitted.at(r_order.at(i)) = true;
  //
  //   // Retrieve index of x_i and x_j
  //   // if (i < 2){
  //   //   int from = INDEX_FROM_KNN(i+1, k);
  //   //   Rcout << "rorder_i: " << r_order.at(i) << std::endl;
  //   //   Rcout << "index from: " << ((int) from) << std::endl;
  //   //   Rcout << "to: " << knn_indices.at( r_order.at(i)) - 1 << std::endl;
  //   //   Rcout << "to_i: " <<  r_order.at(knn_indices.at( r_order.at(i)) - 1) << std::endl;
  //   // }
  //   int x_i = r_order.at(i); //(int) INDEX_FROM_KNN(i+1, k);
  //   int x_j = int(knn_indices.at(x_i) - 1);
  //   // int to = INDEX_TO(r_order.at(i), n), from = INDEX_FROM(r_order.at(i), n, to);
  //
  //
  //   if (admitted.at(x_i) && admitted.at(x_j)) {
  //     if (components.Find(x_i) != components.Find(x_j)){
  //       mst(crow, _) = NumericVector::create(x_i, x_j, *r);
  //       crow++;
  //       components.Union(x_i, x_j);
  //       if (crow == n) break;
  //     }
  //   }
  // }

  return(List::create(_["mst"] = mst, _["admitted"] = admitted));
}

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
    mst(n_edges, _) = NumericVector::create(min_id, c_i, type == 0 ? min * alpha : min);
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
