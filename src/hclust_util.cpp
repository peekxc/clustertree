#include "hclust_util.h"

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


// [[Rcpp::export]]
IntegerMatrix normalizeIndices(const IntegerMatrix& mst){
  IntegerMatrix new_mst = clone(mst);
  IntegerVector from = new_mst.column(0);
  IntegerVector to = new_mst.column(1);
  IntegerVector ids = Rcpp::union_(from, to);
  new_mst(_, 0) = Rcpp::match(from, ids) - 1; // 0-based
  new_mst(_, 1) = Rcpp::match(to, ids) - 1; // 0-based
  return(new_mst);
}

/* mstToHclust
* Given a minimum spanning tree of the columnar form (<from>, <to>, <height>), create a valid hclust object
* using a disjoint-set structure to track components
* Expects 0-based mst indices.
*/
// [[Rcpp::export]]
List mstToHclust(const IntegerMatrix& mst_, const NumericVector& dist){

  // Check to make sure the indices are proper
  const int n = mst_.nrow() + 1;
  const IntegerVector& from_ids = mst_.column(0);
  const IntegerVector& to_ids = mst_.column(1);
  int min_id = std::min((int) min(from_ids), (int) min(to_ids));
  int max_id = std::max((int) max(from_ids), (int) max(to_ids));
  IntegerMatrix mst;
  if (min_id != 0 || max_id != (n - 1)){ mst = normalizeIndices(mst_); }
  else { mst = clone(mst_); }

  // Extract merge heights and associated order of such heights
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
    IntegerVector crow = mst.row(height_order.at(i));
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
    _["order"] = extractOrder(merge), // the check at the beginning of the function should make this safe
    _["labels"] = R_NilValue,
    _["mst"] = mst
  );
  res.attr("class") = "hclust";
  return(res);
}
