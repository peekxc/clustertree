#ifndef METRICS_H
#define METRICS_H

#include <ANN/ANN.h> // ANN typedefs
#include <algorithm>  // std::transform
#include <Rcpp.h> // NumericVector
using namespace Rcpp;

/*  metrics.h
 *  To facilitate the use of arbitrary metrics in the dual tree-associated algorithms, functors that
 *  are able to compute a metric distance are recorded here. These effectively replace the use of 'annDist' and
 *  it's associated macros for the ANN library. Each metric class is derived from an abstract metric class which
 *  creates the interface that every metric functor is required to have, namely:
 *  1) An overloaded function call 'operator()' which recieves as input two ANNpoint's and returns a the distance
 *  between them (of type ANNdist)
 *  2) An overloaded function call 'operator()' which recieves two ANNpoint's, an index 'i', and a distance 'dist', and
 *  returns the incremented distance between 'dist' and the ith component/dimension between the two points. This is required
 *  to facilitate incremental distance calculations and dramatically reduce the overall number of FLOPs in cases where
 *  the data is sparse or high-dimensional. Assumes each component distance is non-negative.
 *  3) A 'finalize' method which transforms a vector of metric distances into their 'final' quantity. This is useful for
 *  metrics where it there exists a monotonic relationship between an intermediate distance value and the final distance
 *  value, and the intermediate distance is used to reduce computation. Examples of this include reserving the
 *  square root calculation until all distances have been computed for euclidean distance, or reserving a unit-conversion
 *  operation for previously normalized coordinate distances (e.g. converting from normalized to WGS84 coordinates). If
 *  not needed, can simply inherit the parent metric class and do nothing.
 */

// Generic metric abstract class to allow polymorphism
struct Metric {
  virtual inline ANNdist operator()(ANNpoint x_i, ANNpoint x_j) = 0;
  virtual inline ANNdist operator()(ANNpoint x_i, ANNpoint x_j, const int i, ANNdist dist) = 0;
  virtual inline void finalize(Rcpp::NumericVector& dist_x) = 0;
  virtual ~Metric() {}
};

// ----------------------------------------------------------
// --- Functors to compute some notion of metric distance ---
// ----------------------------------------------------------

// L_2 norm (squared)
struct L_2 : Metric {
  const unsigned int d; // dimension of the space
  L_2(const int _d) : d(_d) { }
  inline ANNdist operator()(ANNpoint x_i, ANNpoint x_j){
    ANNdist dist = 0;
    for (int i = 0; i < d; ++i) { dist += (x_j[i] - x_i[i]) * (x_j[i] - x_i[i]); }
    return dist;
  }
  // Allow dimension-specific distance calculation for incremental distance updates
  inline ANNdist operator()(ANNpoint x_i, ANNpoint x_j, const int i, ANNdist dist){
    return dist + ((x_j[i] - x_i[i]) * (x_j[i] - x_i[i]));
  }

  // Finalize by taking the square root
  inline void finalize(Rcpp::NumericVector& dist_x){
    std::transform(dist_x.begin(), dist_x.end(), dist_x.begin(), ::sqrt);
  }
};

// L_1 norm
struct L_1 : Metric {
  const unsigned int d; // dimension of the space
  L_1(const int _d) : d(_d) { }
  inline ANNdist operator()(ANNpoint x_i, ANNpoint x_j){
    ANNdist dist = 0;
    for (int i = 0; i < d; ++i) { dist += std::abs(x_j[i] - x_i[i]); }
    return dist;
  }
  // Allow dimension-specific distance calculation for incremental distance updates
  inline ANNdist operator()(ANNpoint x_i, ANNpoint x_j, const int i, ANNdist dist){
    return dist + ((x_j[i] - x_i[i]) >= 0 ? x_j[i] - x_i[i] : x_i[i] - x_j[i]);
  }

  // Do nothing
  inline void finalize(Rcpp::NumericVector& dist_x){}
};

// L_inf norm
struct L_inf : Metric {
  const unsigned int d; // dimension of the space
  L_inf(const int _d) : d(_d) { }
  inline ANNdist operator()(ANNpoint x_i, ANNpoint x_j){
    ANNdist dist = 0;
    for (int i = 0; i < d; ++i) { dist = std::max(dist, std::abs(x_j[i]- x_i[i])); }
    return dist;
  }
  // Allow dimension-specific distance calculation for incremental distance updates
  inline ANNdist operator()(ANNpoint x_i, ANNpoint x_j, const int i, ANNdist dist){
    return std::max(dist, ((x_j[i] - x_i[i]) >= 0 ? x_j[i] - x_i[i] : x_i[i] - x_j[i]));
  }

  // Do nothing
  inline void finalize(Rcpp::NumericVector& dist_x){}
};



struct L_p : Metric {
  const unsigned int d; // dimension of the space
  const unsigned int p; // p norm
  L_p(const unsigned int _d, const unsigned int _p) : d(_d), p(_p) { }
  inline ANNdist operator()(ANNpoint x_i, ANNpoint x_j){
    ANNdist dist = 0;
    for (int i = 0; i < d; ++i) { dist += std::pow(x_j[i] - x_i[i], p); }
    return dist;
  }
  // Allow dimension-specific distance calculation for incremental distance updates
  inline ANNdist operator()(ANNpoint x_i, ANNpoint x_j, const int i, ANNdist dist){
    return dist + std::pow(x_j[i] - x_i[i], p);
  }

  // Raise to the 1/p power
  inline void finalize(NumericVector& dist_x){
    for (NumericVector::iterator dist = dist_x.begin(); dist != dist_x.end(); ++dist){
      *dist = std::pow(*dist, 1/p);
    }
  }
};


// List metric_chooser(std::string metric_name, List options){
//   Rcpp::Function pmatch = Function("pmatch");
//   StringVector metrics = StringVector();
//   metrics.push_back("euclidean");
//   metrics.push_back("manhattan");
//   metrics.push_back("maximum");
//   metrics.push_back("minkowski");
//   const int metric_idx = pmatch(metric_name, metrics);
//   switch(metric_idx){
//
//   }
// }

#endif