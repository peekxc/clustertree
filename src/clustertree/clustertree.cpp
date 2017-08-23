// Includes
#include "clustertree.h"

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
    _["labels"] = R_NilValue,
    _["mst"] = mst
  );
  res.attr("class") = "hclust";
  return(res);
}
