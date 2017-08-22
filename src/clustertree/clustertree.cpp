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
 * Notes: expects 0-based mst indices, and that all edge indices are 0 <= i < n. Is unsafe otherwise!
 */
List mstToHclust(NumericMatrix mst){

  // Extract merge heights and associated order of such heights
  const int n = mst.nrow() + 1;
  NumericVector dist = mst.column(2);
  IntegerVector height_order = order_(dist) - 1;

  // Set up components
  UnionFind components = UnionFind(n);
  IntegerVector comp_index = rep(0, n);
  std::vector<int> assigned = std::vector<int>(n, 0); // 0 to indicate node has not been assigned

  // Go through each row of MST, tracking component indices
  IntegerMatrix merge = IntegerMatrix(n - 1, 2);
  bool from_singleton, to_singleton;
  for (int i = 0; i < n - 1; ++i) {
    if (i % 1000 == 0) Rcpp::checkUserInterrupt();
    NumericVector crow = mst.row(height_order.at(i));
    int from = crow.at(0), to = crow.at(1);
    int from_comp = components.Find(from), to_comp = components.Find(to);
    components.Union(from, to);
    from_singleton = assigned.at(from) == 0, to_singleton = assigned.at(to) == 0;

    // CASE 1: Two singletons a merging
    if (from_singleton && to_singleton) { // Rcout << "1";
      merge.row(i) = IntegerVector::create(-(from + 1), -(to + 1));
      comp_index.at(components.Find(from)) = i;
    }
    // CASE 2: One singleton and one component are merging
    else if (from_singleton || to_singleton) { // Rcout << "2";
      int leaf = from_singleton ? from : to;
      int comp = from_singleton ? to_comp : from_comp;
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
