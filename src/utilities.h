#ifndef CT_UTIL_H
#define CT_UTIL_H

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

//#define NDEBUG 1 // <-- for 'debug' mode, will print out info to R session
#undef NDEBUG // <-- for 'production' mode, will remove IO
#define PROFILING // <-- for 'profile' mode, will give timings of every function being profiled

// Allows indexing lower triangular (dist objects)
#include <math.h>
#define INDEX_TF(N,to,from) (N)*(to) - (to)*(to+1)/2 + (from) - (to) - (1)
#define INDEX_TO(k, n) n - 2 - floor(sqrt(-8*k + 4*n*(n-1)-7)/2.0 - 0.5)
#define INDEX_FROM(k, n, i) k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2
#define INDEX_FROM_KNN(i, k) i % k == 0 ? (i / k) - 1 : int(i / k);

// Simplifies/Shortens iterator definitions
#define VI(x) std::vector<x>::iterator

// Enables to the redirection of messages to the current R session
// UTIL provides a way to input arbitrary code that might be helpful
// only when debugging code
#ifdef NDEBUG
#define R_INFO(x) Rcpp::Rcout << x;
#define R_PRINTF(x, ...) Rprintf(x, __VA_ARGS__);
#define UTIL(x) x;
#else
#define R_INFO(x)
#define R_PRINTF(x, ...)
#define UTIL(x)
#endif

// Trivial profiling technique
#ifdef ANN_PERF
  #include <ANN/ANNperf.h>
  #include <chrono>
  #define BEGIN_PROFILE() { std::chrono::steady_clock::time_point begin_t = std::chrono::steady_clock::now();
  #define END_PROFILE() std::chrono::steady_clock::time_point end_t = std::chrono::steady_clock::now();
  #define REPORT_TIME(x) \
  std::chrono::steady_clock::time_point end_t = std::chrono::steady_clock::now(); \
  Rcout << x << " took: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_t - begin_t).count() << " ms\n"; \
  }
#else
#define BEGIN_PROFILE()
#define END_PROFILE()
#define REPORT_TIME(x)
#endif

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

// TODO: Make compiler-specific options for __restrict__

template <typename T, typename C> bool contains (const T& container, const C& key)
{
  if (std::find(container.begin(), container.end(), key) != container.end()){
    return true;
  } else {
    return false;
  }
}

template<typename T> class CompareIdxByVec { std::vector<T>* _values; public: CompareIdxByVec(std::vector<T>* values) : _values(values) {} public: bool operator() (const size_t& a, const size_t& b) const { return (_values)[a] > (_values)[b]; } };
template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {

  // initialize original index locations
  std::vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(), CompareIdxByVec<T>());

  return idx;
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
  double_edge(){}; // don't initialize anything with default constructor
  double_edge(int from_id, int to_id, double cost) : from(from_id), to(to_id), weight(cost) { }
};

#endif