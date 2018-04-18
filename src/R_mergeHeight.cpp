#include "RcppHeader.h"

// Expects 1-based valid hclust object
// [[Rcpp::export]]
NumericVector mergeHeight(List hcl) {
  IntegerMatrix merge = hcl["merge"];
  NumericVector height = hcl["height"];

  const int n = merge.nrow() + 1;
  NumericVector merge_height = NumericVector(n);

  // Traverse merge
  for (int i = 0; i < n -1; ++i){
    int lm = merge(i, 0), rm = merge(i, 1);
    IntegerVector m = IntegerVector::create(lm, rm);
    if (all(m < 0).is_true()){
      merge_height.at((-m[0]) - 1) = height[i];
      merge_height.at((-m[1]) - 1) = height[i];
    } else if (any(m < 0).is_true()) {
      int leaf = m[0] < 0 ? m[0] : m[1];
      merge_height.at((-leaf) - 1) = height[i];
    }
  }

  // Return ordered merge height
  return merge_height;
}

