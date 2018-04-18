#include "Rcpp.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector test1(const NumericMatrix& x) {
  const int n = x.nrow();
  double dummy_res = 0.0f, x2 = 0.0f;
  std::vector<double> res_store = std::vector<double>();
  for (int i = 0; i < n; ++i){
    for (int j = 0; j < n; ++j){
      dummy_res = (double) sqrt((sqrt(i) + sqrt(j))*(i+j));
      x2 = dummy_res;

      dummy_res = (double) sqrt((sqrt(i) + sqrt(j))*(i+j)) + (i*j)/10.0;
      if ((i % 10) == 0){ res_store.push_back(dummy_res); }

      dummy_res = (double) sqrt((sqrt(i) + sqrt(j))*(i+j)) + + (i*j)/100.0;
      if (i < j && j % 2 == 0){
        std::vector<double>::iterator found_dummy = std::find(res_store.begin(), res_store.end(), dummy_res);
        if (found_dummy != res_store.end()){
          res_store.push_back((double)(i));
          x2 = dummy_res;
        }
      }
    }
  }
  res_store.push_back(x2);
  NumericVector res = wrap(res_store);
  return(res);
}

// [[Rcpp::export]]
NumericVector test2(const NumericMatrix& x) {
  const int n = x.nrow();
  double dummy_res = 0.0f, x2 = 0.0f;
  std::vector<double> res_store = std::vector<double>();


  for (int i = 0; i < n; ++i){
    for (int j = 0; j < n; ++j){
      dummy_res = (double) sqrt((sqrt(i) + sqrt(j))*(i+j));
      x2 = dummy_res;
    }
  }

  for (int i = 0; i < n; ++i){
    for (int j = 0; j < n; ++j){
      dummy_res = (double) sqrt((sqrt(i) + sqrt(j))*(i+j)) + (i*j)/10.0;
      if ((i % 10) == 0){ res_store.push_back(dummy_res); }
    }
  }

  for (int i = 0; i < n; ++i){
    for (int j = 0; j < n; ++j){
      dummy_res = (double) sqrt((sqrt(i) + sqrt(j))*(i+j)) + + (i*j)/100.0;
      if (i < j && j % 2 == 0){
        std::vector<double>::iterator found_dummy = std::find(res_store.begin(), res_store.end(), dummy_res);
        if (found_dummy != res_store.end()){
          res_store.push_back((double)(i));
          x2 = dummy_res;
        }
      }
    }
  }

  res_store.push_back(x2);
  return(wrap(res_store));
}


template<typename Iter, typename Pred, typename Op>
void for_each_if(Iter first, Iter last, Pred p, Op op) {
  while(first != last) {
    if (p(*first)) op(*first);
    ++first;
  }
}

NumericVector test3(const NumericMatrix& x) {
  const int n = x.nrow();
  double dummy_res = 0.0f, x2 = 0.0f;
  std::vector<double> res_store = std::vector<double>();
  for (int i = 0; i < n; ++i){
    for (int j = 0; j < n; ++j){
      dummy_res = (double) sqrt((sqrt(i) + sqrt(j))*(i+j));
      x2 = dummy_res;

      dummy_res = (double) sqrt((sqrt(i) + sqrt(j))*(i+j)) + (i*j)/10.0;
      if ((i % 10) == 0){ res_store.push_back(dummy_res); }

      dummy_res = (double) sqrt((sqrt(i) + sqrt(j))*(i+j)) + + (i*j)/100.0;
      if (i < j && j % 2 == 0){
        std::vector<double>::iterator found_dummy = std::find(res_store.begin(), res_store.end(), dummy_res);
        if (found_dummy != res_store.end()){
          res_store.push_back((double)(i));
          x2 = dummy_res;
        }
      }
    }
  }
  res_store.push_back(x2);
  NumericVector res = wrap(res_store);
  return(res);
}



/*** R
n <- 100
x <- cbind(rnorm(n), rnorm(n))
microbenchmark::microbenchmark({ x1 <- test1(x) }, times = 100L)
microbenchmark::microbenchmark({ x2 <- test2(x) }, times = 100L)
*/
