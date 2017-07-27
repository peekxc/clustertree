#include <Rcpp.h>
using namespace Rcpp;

void int_ptr(std::vector<int>& ids){
  int bkt[6] = {1, 2, 3, 4, 5, 6};
  ids.insert(ids.end(), &bkt[0], &bkt[6]);
}

// [[Rcpp::export]]
IntegerVector test_ptr() {
  std::vector<int> test_vec = std::vector<int>();
  int_ptr(test_vec);
  return(wrap(test_vec));
}

/*** R
test_ptr()
  */
