#ifndef CT_UTIL_H
#define CT_UTIL_H

#include <Rcpp.h>
using namespace Rcpp;

// Allows indexing lower triangular (dist objects)
#include <math.h>
#define INDEX_TF(N,to,from) (N)*(to) - (to)*(to+1)/2 + (from) - (to) - (1)
#define INDEX_TO(k, n) n - 2 - floor(sqrt(-8*k + 4*n*(n-1)-7)/2.0 - 0.5)
#define INDEX_FROM(k, n, i) k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2
#define INDEX_FROM_KNN(i, k) i % k == 0 ? (i / k) - 1 : int(i / k);


// std::to_string is apparently a c++11 only thing that crashes appveyor, so using ostringstream it is!
namespace patch
{
template < typename T > std::string to_string( const T& n )
{
  std::ostringstream stm ;
  stm << n ;
  return stm.str() ;
}
}

template <typename T, typename C> bool contains (const T& container, const C& key)
{
  if (std::find(container.begin(), container.end(), key) != container.end()){
    return true;
  } else {
    return false;
  }
}

inline IntegerVector lowerTri(IntegerMatrix m){
  int n = m.nrow();
  IntegerVector lower_tri = IntegerVector(n * (n - 1) / 2);
  for (int i = 0, c = 0; i < n; ++i){
    for (int j = i + 1; j < n; ++j){
      if (i < j) lower_tri[c++] = m(i, j);
    }
  }
  return(lower_tri);
}

inline IntegerVector which_cpp(LogicalVector x, bool value) {
  int nx = x.size();
  std::vector<int> y;
  y.reserve(nx);
  for(int i = 0; i < nx; ++i) { if (x[i] == value) y.push_back(i); }
  return wrap(y);
}


inline IntegerVector which_cpp( NumericVector x, double value) {
  int nx = x.size();
  std::vector<int> y;
  y.reserve(nx);
  for(int i = 0; i < nx; ++i) { if (x[i] == value) y.push_back(i); }
  return wrap(y);
}

inline IntegerVector which_cpp( IntegerVector x, int value) {
  int nx = x.size();
  std::vector<int> y;
  y.reserve(nx);
  for(int i = 0; i < nx; ++i) { if (x[i] == value) y.push_back(i); }
  return wrap(y);
}

inline IntegerVector which_geq( IntegerVector x, int value) {
  int nx = x.size();
  std::vector<int> y;
  y.reserve(nx);
  for(int i = 0; i < nx; ++i) { if (x[i] >= value) y.push_back(i); }
  return wrap(y);
}

inline NumericVector combine(const NumericVector& t1, const NumericVector& t2){
  std::size_t n = t1.size() + t2.size();
  NumericVector output = Rcpp::no_init(n);
  std::copy(t1.begin(), t1.end(), output.begin());
  std::copy(t2.begin(), t2.end(), output.begin()+t1.size());
  return output;
}

inline IntegerVector combine(const IntegerVector& t1, const IntegerVector& t2){
  std::size_t n = t1.size() + t2.size();
  IntegerVector output = Rcpp::no_init(n);
  std::copy(t1.begin(), t1.end(), output.begin());
  std::copy(t2.begin(), t2.end(), output.begin()+t1.size());
  return output;
}

// Faster version of above combine function, assuming you can precompute and store
// the containers needing to be concatenated
inline IntegerVector concat_int (List const& container){
  int total_length = 0;
  for (List::const_iterator it = container.begin(); it != container.end(); ++it){
      total_length += as<IntegerVector>(*it).size();
  }
  int pos = 0;
  IntegerVector output = Rcpp::no_init(total_length);
  for (List::const_iterator it = container.begin(); it != container.end(); ++it){
    IntegerVector vec = as<IntegerVector>(*it);
    std::copy(vec.begin(), vec.end(), output.begin() + pos);
    pos += vec.size();
  }
  return(output);
}

// Based on (but extended) http://stackoverflow.com/questions/21609934/ordering-permutation-in-rcpp-i-e-baseorder
inline IntegerVector order_(NumericVector x) {
  if (is_true(any(duplicated(x)))) {
    Function order("order");
    return(as<IntegerVector>(order(x)));
  }
  NumericVector sorted = clone(x).sort();
  return match(sorted, x);
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

struct double_edge
{
  unsigned int from, to;
  double weight;
  double_edge(int from_id, int to_id, double cost) : from(from_id), to(to_id), weight(cost) { }
};

#endif