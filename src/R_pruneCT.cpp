#include <Rcpp.h>
using namespace Rcpp;

#include "union_find.h"

// [[Rcpp::export]]
NumericVector pruneCT(List C_n, NumericVector prune_heights, IntegerVector valid_idx) {

  IntegerMatrix merge = C_n["merge"];
  NumericVector height = C_n["height"];

  // Create
  const int n = height.length() + 1;
  UnionFind old_comp = UnionFind(n), new_comp = UnionFind(n);

  // Recreate component structure, updating
  IntegerMatrix pruned_merge = IntegerMatrix(n - 1, 2), old_merge = IntegerMatrix(n - 1, 2);
  IntegerVector comp_index = rep(0, n);
  std::vector<int> assigned = std::vector<int>(n, 0); // 0 to indicate node has not been assigned


  // The current index to the height r' > r
  IntegerVector::iterator prune_idx = valid_idx.begin();
  for (int i = 0, n_i = 0; i < n - 1; ++i) {
    if (i % 1000 == 0) Rcpp::checkUserInterrupt();
    while(n_i < *prune_idx){

    }
    if (prune_idx != valid_idx.end() && i == *prune_idx){
      // Pruning step: Connect any two components of C_n(r) that belong to the
      // same connected component in C_n(r(max (λ̃r, 0))).
      for (int cc = 0; cc < n; ++cc){ new_comp.Union(new_comp.Find(cc), old_comp.Find(cc)); }
      prune_idx++;
    }

    IntegerVector crow = merge.row(i);
    int from = crow.at(0), to = crow.at(1);
    if (from < 0 && to < 0) {
      old_merge.row(i) = IntegerVector::create(from, to);
      //pruned_merge.row(i) = IntegerVector::create(from, to);
      // comp_index.at(i)
    }
    else if (from < 0 || to < 0) {
      int leaf = from < 0 ? from: to;
      int comp = from < 0 ? to : from;
      merge.row(i) = IntegerVector::create(-(leaf), comp_index.at(comp)+1);
      //comp_index.at(components.Find(leaf)) = i;
    }
  }


  //return x * 2;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
timesTwo(42)
*/
